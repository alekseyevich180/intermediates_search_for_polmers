import argparse
import time
import weakref
from pathlib import Path
from typing import IO, Any, Union

import numpy as np
from ase import Atoms, units
from ase.calculators.calculator import PropertyNotImplementedError
from ase.eos import EquationOfState
from ase.io import read, write
from ase.md import MDLogger
from ase.md.npt import NPT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.parallel import world
from ase.utils import IOContext
from fairchem.core import pretrained_mlip, FAIRChemCalculator
import os
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

# print(f"pfp_api_client: {pfp_api_client.__version__}")

import io
import json
import pandas as pd
from ase.constraints import ExpCellFilter, StrainFilter
from ase.io.jsonio import write_json, read_json
from ase.optimize import LBFGS, FIRE
from IPython.display import Image
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import optuna
from ase.visualize import view
from ase.neb import NEB
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from ase.geometry import get_distances
from ase.constraints import FixAtoms
from ase.optimize import BFGSLineSearch
# from scipy.spatial.distance import pdist, squareform  # 暂时注释掉
# from sklearn.cluster import DBSCAN  # 暂时注释掉，避免导入错误

from pfcc_extras.visualize.view import view_ngl
from pfcc_extras.visualize.ase import view_ase_atoms


class FunctionalGroupAnalyzer:
    """分析分子中的官能团和周围的C原子"""
    
    def __init__(self, atoms, calculator):
        self.atoms = atoms
        self.calculator = calculator
        self.functional_groups = []
        self.carbon_atoms = []
        self.bond_network = {}
        
    def identify_functional_groups(self):
        """识别分子中的官能团"""
        symbols = self.atoms.get_chemical_symbols()
        positions = self.atoms.get_positions()
        
        # 定义官能团模式
        functional_patterns = {
            'ketone': ['C', 'O'],  # 酮基
            'alcohol': ['O', 'H'],  # 羟基
            'carboxylic_acid': ['C', 'O', 'O'],  # 羧基
            'aldehyde': ['C', 'O'],  # 醛基
            'ester': ['C', 'O', 'C'],  # 酯基
            'amine': ['N', 'H'],  # 氨基
            'amide': ['C', 'O', 'N'],  # 酰胺基
        }
        
        # 查找官能团
        for group_name, pattern in functional_patterns.items():
            group_atoms = self._find_pattern(pattern, symbols, positions)
            if group_atoms:
                self.functional_groups.append({
                    'name': group_name,
                    'atoms': group_atoms,
                    'center': self._calculate_group_center(group_atoms, positions)
                })
        
        return self.functional_groups
    
    def _find_pattern(self, pattern, symbols, positions):
        """查找特定的原子模式"""
        # 简化的模式匹配，实际应用中可能需要更复杂的算法
        pattern_atoms = []
        
        if pattern == ['C', 'O']:  # 酮基或醛基
            for i, symbol in enumerate(symbols):
                if symbol == 'C':
                    for j, other_symbol in enumerate(symbols):
                        if other_symbol == 'O' and i != j:
                            distance = np.linalg.norm(positions[i] - positions[j])
                            if distance < 1.5:  # C=O键长
                                pattern_atoms = [i, j]
                                break
                    if pattern_atoms:
                        break
        
        return pattern_atoms
    
    def _calculate_group_center(self, atom_indices, positions):
        """计算官能团的几何中心"""
        if not atom_indices:
            return None
        return np.mean([positions[i] for i in atom_indices], axis=0)
    
    def find_carbons_near_functional_groups(self, max_distance=3.0):
        """找到官能团周围的C原子"""
        if not self.functional_groups:
            self.identify_functional_groups()
        
        symbols = self.atoms.get_chemical_symbols()
        positions = self.atoms.get_positions()
        
        nearby_carbons = []
        
        for group in self.functional_groups:
            group_center = group['center']
            if group_center is None:
                continue
            
            for i, symbol in enumerate(symbols):
                if symbol == 'C' and i not in group['atoms']:  # 排除官能团中的C
                    distance = np.linalg.norm(positions[i] - group_center)
                    if distance <= max_distance:
                        nearby_carbons.append({
                            'index': i,
                            'distance_to_group': distance,
                            'functional_group': group['name'],
                            'position': positions[i].copy()
                        })
        
        # 按距离排序
        nearby_carbons.sort(key=lambda x: x['distance_to_group'])
        self.carbon_atoms = nearby_carbons
        
        return nearby_carbons
    
    def build_bond_network(self):
        """构建键网络，用于约束优化"""
        positions = self.atoms.get_positions()
        symbols = self.atoms.get_chemical_symbols()
        
        bond_network = {}
        
        for i in range(len(self.atoms)):
            bond_network[i] = []
            for j in range(i+1, len(self.atoms)):
                distance = np.linalg.norm(positions[i] - positions[j])
                # 判断是否为键（简化判断）
                if self._is_bond(i, j, distance, symbols):
                    bond_network[i].append(j)
                    bond_network[j].append(i)
        
        self.bond_network = bond_network
        return bond_network
    
    def _is_bond(self, i, j, distance, symbols):
        """判断两个原子是否成键"""
        # 简化的键判断逻辑
        bond_lengths = {
            ('C', 'C'): 1.6,
            ('C', 'O'): 1.5,
            ('C', 'N'): 1.5,
            ('C', 'H'): 1.2,
            ('O', 'H'): 1.0,
            ('N', 'H'): 1.1,
        }
        
        pair = tuple(sorted([symbols[i], symbols[j]]))
        max_bond_length = bond_lengths.get(pair, 2.0)
        
        return distance < max_bond_length


class MDStabilitySearcher:
    """通过MD模拟寻找最稳定的结构"""
    
    def __init__(self, atoms, calculator, temperature=300, steps=100000, sample_interval=100):
        self.atoms = atoms
        self.calculator = calculator
        self.temperature = temperature
        self.steps = steps
        self.sample_interval = sample_interval
        self.sampled_structures = []
        self.energies = []
        
    def run_md_sampling(self):
        """运行MD模拟并采样结构"""
        print(f"🧪 Running MD sampling: {self.steps} steps at {self.temperature}K")
        
        # 设置计算器
        self.atoms.calc = self.calculator
        
        # 设置初始速度
        MaxwellBoltzmannDistribution(self.atoms, temperature_K=self.temperature, force_temp=True)
        Stationary(self.atoms)
        
        # 创建MD动力学
        dyn = NPT(
            self.atoms,
            0.5 * units.fs,
            temperature_K=self.temperature,
            externalstress=0,
            ttime=100 * units.fs,
            pfactor=None,
            trajectory="md_sampling.xyz",
            loginterval=self.sample_interval,
        )
        
        # 采样回调函数
        def sample_structure():
            step = dyn.get_number_of_steps()
            if step % self.sample_interval == 0:
                energy = self.atoms.get_potential_energy()
                structure = self.atoms.copy()
                
                self.sampled_structures.append({
                    'step': step,
                    'atoms': structure,
                    'energy': energy,
                    'positions': structure.get_positions().copy()
                })
                self.energies.append(energy)
                
                if step % (self.sample_interval * 10) == 0:
                    print(f"   Sampled at step {step}, energy: {energy:.4f} eV")
        
        # 附加采样
        dyn.attach(sample_structure, interval=self.sample_interval)
        
        # 运行MD
        dyn.run(self.steps)
        
        print(f"   Sampling complete: {len(self.sampled_structures)} structures sampled")
        return self.sampled_structures
    
    def find_most_stable_structure(self):
        """找到最稳定的结构"""
        if not self.sampled_structures:
            self.run_md_sampling()
        
        # 找到能量最低的结构
        min_energy_idx = np.argmin(self.energies)
        most_stable = self.sampled_structures[min_energy_idx]
        
        print(f"🎯 Most stable structure found at step {most_stable['step']}")
        print(f"   Energy: {most_stable['energy']:.4f} eV")
        
        return most_stable


class ProgressiveCarbonOptimizer:
    """从近到远逐步优化C原子的位置"""
    
    def __init__(self, atoms, calculator, carbon_atoms, functional_groups):
        self.atoms = atoms
        self.calculator = calculator
        self.carbon_atoms = carbon_atoms
        self.functional_groups = functional_groups
        self.optimization_history = []
        
    def optimize_carbons_progressively(self, n_trials=100):
        """逐步优化C原子位置"""
        print(f"🔧 Starting progressive carbon optimization...")
        print(f"   {len(self.carbon_atoms)} carbon atoms to optimize")
        
        # 创建Optuna研究
        study = optuna.create_study(direction='minimize')
        
        # 设置用户属性
        study.set_user_attr("atoms", atoms_to_json(self.atoms))
        study.set_user_attr("calculator", "ASECalculator")
        study.set_user_attr("carbon_atoms", self.carbon_atoms)
        
        # 逐步优化每个C原子
        for i, carbon_info in enumerate(self.carbon_atoms):
            print(f"\n🎯 Optimizing carbon atom {i+1}/{len(self.carbon_atoms)}")
            print(f"   Distance to functional group: {carbon_info['distance_to_group']:.2f} Å")
            
            # 为当前C原子创建优化目标函数
            def objective(trial):
                return self._optimize_single_carbon(trial, carbon_info, i)
            
            # 运行优化
            study.optimize(objective, n_trials=n_trials)
            
            # 保存最佳结果
            best_params = study.best_params
            best_energy = study.best_value
            
            print(f"   Best energy: {best_energy:.4f} eV")
            print(f"   Best parameters: {best_params}")
            
            # 更新原子结构
            self._apply_optimization_result(best_params, carbon_info)
            
            # 记录历史
            self.optimization_history.append({
                'carbon_index': i,
                'carbon_info': carbon_info,
                'best_energy': best_energy,
                'best_params': best_params
            })
        
        return self.optimization_history
    
    def _optimize_single_carbon(self, trial, carbon_info, carbon_index):
        """优化单个C原子的位置"""
        # 创建原子副本
        atoms_copy = self.atoms.copy()
        atoms_copy.calc = self.calculator
        
        # 固定其他原子（除了当前优化的C原子）
        carbon_idx = carbon_info['index']
        constraints = []
        
        for i in range(len(atoms_copy)):
            if i != carbon_idx:
                constraints.append(FixAtoms(indices=[i]))
        
        atoms_copy.set_constraint(constraints)
        
        # 获取当前C原子的位置
        current_pos = atoms_copy.get_positions()[carbon_idx]
        
        # 定义优化参数（位置和旋转）
        # 位置偏移
        dx = trial.suggest_float(f"dx_{carbon_index}", -1.0, 1.0)
        dy = trial.suggest_float(f"dy_{carbon_index}", -1.0, 1.0)
        dz = trial.suggest_float(f"dz_{carbon_index}", -0.5, 0.5)
        
        # 旋转角度
        phi = trial.suggest_float(f"phi_{carbon_index}", 0, 2*np.pi)
        theta = trial.suggest_float(f"theta_{carbon_index}", 0, np.pi)
        psi = trial.suggest_float(f"psi_{carbon_index}", 0, 2*np.pi)
        
        # 应用位置偏移
        new_pos = current_pos + np.array([dx, dy, dz])
        atoms_copy.positions[carbon_idx] = new_pos
        
        # 检查键约束
        if not self._check_bond_constraints(atoms_copy, carbon_idx):
            return float('inf')  # 违反键约束，返回无穷大能量
        
        # 优化结构
        try:
            opt = LBFGS(atoms_copy)
            opt.run(fmax=0.01)
            
            energy = atoms_copy.get_potential_energy()
            
            # 再次检查键约束
            if not self._check_bond_constraints(atoms_copy, carbon_idx):
                return float('inf')
            
            return energy
            
        except Exception as e:
            print(f"   Optimization failed: {e}")
            return float('inf')
    
    def _check_bond_constraints(self, atoms, carbon_idx):
        """检查键约束是否满足"""
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        # 检查当前C原子与相邻原子的键长
        for i in range(len(atoms)):
            if i != carbon_idx:
                distance = np.linalg.norm(positions[carbon_idx] - positions[i])
                symbols_pair = tuple(sorted([symbols[carbon_idx], symbols[i]]))
                
                # 定义键长范围
                bond_ranges = {
                    ('C', 'C'): (1.2, 1.8),
                    ('C', 'O'): (1.1, 1.7),
                    ('C', 'N'): (1.1, 1.7),
                    ('C', 'H'): (0.9, 1.3),
                }
                
                if symbols_pair in bond_ranges:
                    min_bond, max_bond = bond_ranges[symbols_pair]
                    if distance < min_bond or distance > max_bond:
                        return False  # 键长超出合理范围
        
        return True
    
    def _apply_optimization_result(self, params, carbon_info):
        """应用优化结果到原子结构"""
        carbon_idx = carbon_info['index']
        current_pos = self.atoms.get_positions()[carbon_idx]
        
        # 提取参数
        dx = params.get(f"dx_{self.carbon_atoms.index(carbon_info)}", 0)
        dy = params.get(f"dy_{self.carbon_atoms.index(carbon_info)}", 0)
        dz = params.get(f"dz_{self.carbon_atoms.index(carbon_info)}", 0)
        
        # 应用位置变化
        new_pos = current_pos + np.array([dx, dy, dz])
        self.atoms.positions[carbon_idx] = new_pos
        
        # 重新优化整个结构
        self.atoms.calc = self.calculator
        opt = LBFGS(self.atoms)
        opt.run(fmax=0.01)


class ReactionDetector:
    """检测和分析氧化反应中的中间状态"""
    
    def __init__(self, atoms, calculator, bond_threshold=2.0, save_interval=100):
        self.atoms = atoms
        self.calculator = calculator
        self.bond_threshold = bond_threshold  # 键长阈值 (Å)
        self.save_interval = save_interval
        self.intermediates = []
        self.bond_changes = []
        self.energy_history = []
        self.reaction_steps = []
        
    def analyze_bonds(self, atoms):
        """分析分子中的化学键"""
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        bonds = []
        
        for i in range(len(atoms)):
            for j in range(i+1, len(atoms)):
                distance = np.linalg.norm(positions[i] - positions[j])
                if distance < self.bond_threshold:
                    bond_type = f"{symbols[i]}-{symbols[j]}"
                    bonds.append({
                        'atoms': (i, j),
                        'distance': distance,
                        'type': bond_type,
                        'symbols': (symbols[i], symbols[j])
                    })
        return bonds
    
    def detect_bond_changes(self, bonds_prev, bonds_curr):
        """检测键的变化"""
        changes = []
        
        # 检测新形成的键
        for bond in bonds_curr:
            found = False
            for prev_bond in bonds_prev:
                if set(bond['atoms']) == set(prev_bond['atoms']):
                    found = True
                    # 检查键长变化
                    if abs(bond['distance'] - prev_bond['distance']) > 0.1:
                        changes.append({
                            'type': 'stretch',
                            'atoms': bond['atoms'],
                            'old_distance': prev_bond['distance'],
                            'new_distance': bond['distance'],
                            'bond_type': bond['type']
                        })
                    break
            if not found:
                changes.append({
                    'type': 'formation',
                    'atoms': bond['atoms'],
                    'distance': bond['distance'],
                    'bond_type': bond['type']
                })
        
        # 检测断裂的键
        for bond in bonds_prev:
            found = False
            for curr_bond in bonds_curr:
                if set(bond['atoms']) == set(curr_bond['atoms']):
                    found = True
                    break
            if not found:
                changes.append({
                    'type': 'breaking',
                    'atoms': bond['atoms'],
                    'old_distance': bond['distance'],
                    'bond_type': bond['type']
                })
        
        return changes
    
    def is_reaction_intermediate(self, atoms, step, energy):
        """判断是否为反应中间体"""
        bonds = self.analyze_bonds(atoms)
        
        # 如果有键的变化，可能是中间体
        if len(self.bond_changes) > 0:
            bond_changes = self.detect_bond_changes(
                self.intermediates[-1]['bonds'] if self.intermediates else [],
                bonds
            )
            
            if bond_changes:
                intermediate = {
                    'step': step,
                    'atoms': atoms.copy(),
                    'energy': energy,
                    'bonds': bonds,
                    'bond_changes': bond_changes,
                    'coordinates': atoms.get_positions().copy(),
                    'symbols': atoms.get_chemical_symbols()
                }
                
                self.intermediates.append(intermediate)
                self.bond_changes.extend(bond_changes)
                self.reaction_steps.append(step)
                
                return True, intermediate
        
        return False, None
    
    def save_intermediate(self, intermediate, output_dir="intermediates"):
        """保存中间体结构"""
        os.makedirs(output_dir, exist_ok=True)
        
        step = intermediate['step']
        atoms = intermediate['atoms']
        
        # 保存结构文件
        write(f"{output_dir}/intermediate_step_{step:06d}.cif", atoms)
        
        # 保存详细信息
        info = {
            'step': step,
            'energy': intermediate['energy'],
            'bond_changes': intermediate['bond_changes'],
            'bonds': intermediate['bonds']
        }
        
        with open(f"{output_dir}/intermediate_step_{step:06d}.json", 'w') as f:
            json.dump(info, f, indent=2, default=str)
    
    def analyze_reaction_path(self):
        """分析反应路径"""
        if len(self.intermediates) < 2:
            return None
        
        path_info = {
            'total_intermediates': len(self.intermediates),
            'reaction_steps': self.reaction_steps,
            'energy_profile': [int['energy'] for int in self.intermediates],
            'bond_changes_summary': {}
        }
        
        # 统计键变化
        for change in self.bond_changes:
            bond_type = change['bond_type']
            if bond_type not in path_info['bond_changes_summary']:
                path_info['bond_changes_summary'][bond_type] = []
            path_info['bond_changes_summary'][bond_type].append(change['type'])
        
        return path_info


class IntermediateTracker:
    """跟踪和保存反应中间体"""
    
    def __init__(self, output_dir="reaction_intermediates"):
        self.output_dir = output_dir
        self.intermediates = []
        self.energy_threshold = 0.1  # eV
        self.geometry_threshold = 0.5  # Å
        
    def should_save_intermediate(self, atoms, energy):
        """判断是否应该保存这个中间体"""
        if not self.intermediates:
            return True
        
        # 检查能量差异
        last_energy = self.intermediates[-1]['energy']
        if abs(energy - last_energy) < self.energy_threshold:
            return False
        
        # 检查几何差异
        last_positions = self.intermediates[-1]['positions']
        current_positions = atoms.get_positions()
        
        # 计算RMSD
        rmsd = np.sqrt(np.mean((current_positions - last_positions)**2))
        if rmsd < self.geometry_threshold:
            return False
        
        return True
    
    def add_intermediate(self, atoms, energy, step, bonds=None, forces=None):
        """添加中间体"""
        if self.should_save_intermediate(atoms, energy):
            intermediate = {
                'step': step,
                'atoms': atoms.copy(),
                'energy': energy,
                'positions': atoms.get_positions().copy(),
                'symbols': atoms.get_chemical_symbols(),
                'bonds': bonds,
                'forces': forces
            }
            
            self.intermediates.append(intermediate)
            self.save_intermediate(intermediate)
    
    def save_intermediate(self, intermediate):
        """保存中间体"""
        step = intermediate['step']
        
        # 保存结构
        write(f"{self.output_dir}/step_{step:06d}.cif", intermediate['atoms'])
        
        # 保存能量和几何信息
        info = {
            'step': step,
            'energy': intermediate['energy'],
            'positions': intermediate['positions'].tolist(),
            'symbols': intermediate['symbols']
        }
        
        if intermediate['bonds']:
            info['bonds'] = intermediate['bonds']
        
        with open(f"{self.output_dir}/step_{step:06d}.json", 'w') as f:
            json.dump(info, f, indent=2)


class MDLogger(IOContext):
    """Class for logging molecular dynamics simulations."""

    def __init__(
        self,
        dyn: Any,  # not fully annotated so far to avoid a circular import
        atoms: Atoms,
        logfile: Union[IO, str],
        stress: bool = False,
        mode: str = "a",
        comm=world,
    ):
        """
        Args:
            dyn (Any): The dynamics.  Only a weak reference is kept.
            atoms (Atoms): The atoms.
            logfile (Union[IO, str]): File name or open file, "-" meaning standard output.
            stress (bool, optional): Include stress in log.
            mode (str, optional): How the file is opened if logfile is a filename.
        """
        self.dyn = weakref.proxy(dyn) if hasattr(dyn, "get_time") else None
        self.atoms = atoms
        global_natoms = atoms.get_global_number_of_atoms()
        self.logfile = self.openfile(file=logfile, mode=mode, comm=comm)
        self.stress = stress
        self.hdr = "%-9s %-9s %7s %12s %12s %12s  %12s" % (
            "Step",
            "Time[ps]",
            "T[K]",
            "Epot[eV]",
            "Ekin[eV]",
            "Etot[eV]",
            "Density[g/cm3]",
        )
        # Choose a sensible number of decimals
        if global_natoms <= 100:
            digits = 4
        elif global_natoms <= 1000:
            digits = 3
        elif global_natoms <= 10000:
            digits = 2
        else:
            digits = 1
        self.fmt = "%-10d %-10.4f %6.1f" + 4 * ("%%12.%df " % (digits))
        if self.stress:
            self.hdr += "   ----------------------- stress [GPa] ------------------------"
            self.fmt += 6 * " %10.3f"
        self.fmt += "\n"
        self.logfile.write(self.hdr + "\n")

    def __del__(self):
        self.close()

    def __call__(self):
        if self.dyn:
            t = self.dyn.get_time() / (1000 * units.fs)
            temp = self.atoms.get_temperature()
            epot = self.atoms.get_potential_energy()
            ekin = self.atoms.get_kinetic_energy()
            density = (sum(self.atoms.get_masses()) / units.mol) / (self.atoms.get_volume() * 1e-24)
            dat = (self.dyn.nsteps, t, temp, epot, ekin, epot + ekin, density)
            if self.stress:
                dat += tuple(self.atoms.get_stress(include_ideal_gas=True) / units.GPa)
            self.logfile.write(self.fmt % dat)
            self.logfile.flush()


def shape_upper_triangular_cell(atoms: Atoms) -> Atoms:
    """Transform to upper-triangular cell.

    Args:
        atoms (Atoms): Atoms objects

    Returns:
        atoms (Atoms): Atoms objects whose cell is shaped to upper-triangular cell
    """
    if not NPT._isuppertriangular(atoms.get_cell()):
        a, b, c, alpha, beta, gamma = atoms.cell.cellpar()
        angles = np.radians((alpha, beta, gamma))
        sin_a, sin_b, _sin_g = np.sin(angles)
        cos_a, cos_b, cos_g = np.cos(angles)
        cos_p = (cos_g - cos_a * cos_b) / (sin_a * sin_b)
        cos_p = np.clip(cos_p, -1, 1)
        sin_p = (1 - cos_p**2) ** 0.5

        new_basis = [
            (a * sin_b * sin_p, a * sin_b * cos_p, a * cos_b),
            (0, b * sin_a, b * cos_a),
            (0, 0, c),
        ]
        atoms.set_cell(new_basis, scale_atoms=True)
    return atoms


def calculate_eos(atoms: Atoms, linspace_step: int) -> tuple[np.ndarray, np.ndarray]:
    """Calculate energy and volume for given atoms.

    Args:
        atoms (Atoms): ASE atoms object
        linspace_step (int, optional): define interval of volumes. Defaults to 20.

    Returns:
        volumes (np.ndarray): volumes of system (range: 0.95 to 1.05)
        energies (np.ndarray): inferenced energies
    """
    volumes = []
    energies = []
    base_cell = atoms.get_cell()
    # DFT(scf) or NNP inference with different sizes of cell
    for x in np.linspace(0.95, 1.05, linspace_step):
        atoms.set_cell(base_cell * x, scale_atoms=True)
        volume = atoms.get_volume()
        energy = atoms.get_potential_energy() / len(atoms)
        volumes.append(volume)
        energies.append(energy)

    return np.array(volumes), np.array(energies)

def write_poscar(atoms, step, structure_dir="structures"):
        """保存POSCAR文件"""
        os.makedirs(structure_dir, exist_ok=True)
        fname = os.path.join(structure_dir, f"POSCAR_step{step:06d}.vasp")
        write(fname, atoms, format="vasp")


def create_reaction_monitor(atoms, calculator, output_dir="reaction_analysis"):
    """创建反应监控器"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 创建反应检测器和中间体跟踪器
    detector = ReactionDetector(atoms, calculator, bond_threshold=2.0)
    tracker = IntermediateTracker(output_dir)
    
    return detector, tracker


def monitor_reaction_step(atoms, step, detector, tracker):
    """监控每一步的反应状态"""
    energy = atoms.get_potential_energy()
    
    # 分析键的变化
    bonds = detector.analyze_bonds(atoms)
    
    # 检测是否为反应中间体
    is_intermediate, intermediate = detector.is_reaction_intermediate(atoms, step, energy)
    
    if is_intermediate:
        print(f"🔬 Reaction intermediate detected at step {step}")
        print(f"   Energy: {energy:.4f} eV")
        print(f"   Bond changes: {len(intermediate['bond_changes'])}")
        
        # 保存中间体
        detector.save_intermediate(intermediate, "intermediates")
    
    # 添加到跟踪器
    tracker.add_intermediate(atoms, energy, step, bonds)
    
    return is_intermediate, intermediate


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--atoms_path", type=Path, required=True)
    parser.add_argument("--model_path", type=Path, default=None)
    parser.add_argument("--package", type=str, default=None)
    parser.add_argument("--out_traj_path", type=Path, required=True)
    parser.add_argument("--temperature", type=float, default=300)
    parser.add_argument("--timestep", type=float, default=0.5)
    parser.add_argument("--run_steps", type=int, default=10000)
    parser.add_argument("--traj_interval", type=int, default=100)
    parser.add_argument("--ensemble", type=str, default="nvt", choices=["nvt", "npt"])
    parser.add_argument("--taut", type=int, default=100)
    parser.add_argument("--pressure", type=float, default=1.0)
    parser.add_argument("--taup", type=int, default=1000)
    parser.add_argument("--ensemble_model_paths", type=Path, nargs="*")
    parser.add_argument("--include_d3", action="store_true")
    parser.add_argument("--uma_model", type=str, default="uma-m-1p1",
                    help="UMA model name, e.g., uma-s-1p1 / uma-m-1p1")
    parser.add_argument("--device", type=str, default="cuda",
                    choices=["cuda", "cpu"], help="inference device")
    parser.add_argument("--save_interval", type=int, default=1000)
    args = parser.parse_args()

    print(args)
    args.out_traj_path.parent.mkdir(exist_ok=True, parents=True)
    atoms = read(args.atoms_path)
    atoms = shape_upper_triangular_cell(atoms)

    ### SET CALCULATOR ###
    
    estimator = Estimator(model_version="v8.0.0",calc_mode=EstimatorCalcMode.CRYSTAL_U0_PLUS_D3)
    calculator = ASECalculator(estimator)
    ######################

    if args.include_d3:
        ### DFTD3 implemented in ASE and https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/
        from ase.calculators.dftd3 import DFTD3

        calc = DFTD3(dft=calc)

        ### TorchDFTD3 ###
        # from ase.calculators.mixing import SumCalculator
        # from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
        # d3 = TorchDFTD3Calculator(atoms=atoms, device="cpu")
        # calc = SumCalculator([calc, d3])
    atoms.calc = calculator

    # set momenta
    MaxwellBoltzmannDistribution(atoms, temperature_K=args.temperature, force_temp=True)
    Stationary(atoms)

    # 创建反应监控器
    detector, tracker = create_reaction_monitor(atoms, calculator, "reaction_analysis")
    print("🔬 Reaction monitoring initialized")

    if args.ensemble == "nvt":
        dyn = NPT(
            atoms,
            args.timestep * units.fs,
            temperature_K=args.temperature,
            externalstress=0,
            ttime=args.taut * units.fs,
            pfactor=None,
            trajectory=str(args.out_traj_path),
            loginterval=args.traj_interval,
        )
        logger = MDLogger(dyn, atoms, logfile="-", mode="w", stress=False)
        dyn.attach(lambda: write_poscar(atoms, dyn.get_number_of_steps()), interval=args.save_interval)    
    elif args.ensemble == "npt":
        if args.ensemble_model_paths:
            calc.calculate_model_devi = False
        # Check whether calc can compute stress property
        _atoms = atoms.copy()
        _atoms.calc = calc
        try:
            _atoms.get_stress()
        except PropertyNotImplementedError:
            print("Calculator cannot compute stress property")
            return

        # compute bulk modulus
        print("Computing bulk modulus...")
        v, e = calculate_eos(_atoms, 100)
        eos = EquationOfState(v, e, eos="murnaghan")
        _, _, B = eos.fit()
        bulk_modulus_GPa = B / units.kJ * 1e24
        print(f"Bulk Modulus: {bulk_modulus_GPa} GPa")

        if args.ensemble_model_paths:
            calc.calculate_model_devi = True
        dyn = NPT(
            atoms,
            args.timestep * units.fs,
            temperature_K=args.temperature,
            externalstress=args.pressure * units.bar,
            ttime=args.taut * units.fs,
            pfactor=(args.taup * units.fs) ** 2 * bulk_modulus_GPa * units.GPa,
            mask=np.identity(3),
            trajectory=str(args.out_traj_path),
            loginterval=args.traj_interval,
        )
        logger = MDLogger(dyn, atoms, logfile="-", mode="a", stress=True)
    else:
        print(f"Unsupported ensemble: {args.ensemble}")
        return

    dyn.attach(logger, interval=args.traj_interval)

    # 创建反应监控回调函数
    def reaction_monitor_callback():
        step = dyn.get_number_of_steps()
        is_intermediate, intermediate = monitor_reaction_step(atoms, step, detector, tracker)
        if is_intermediate:
            print(f"   Bond changes: {[change['type'] + ' ' + change['bond_type'] for change in intermediate['bond_changes']]}")

    # 附加反应监控
    dyn.attach(reaction_monitor_callback, interval=args.traj_interval)

    # run md simulation
    time_start = time.perf_counter()
    print("Starting MD simulation with reaction monitoring...")
    dyn.run(args.run_steps)
    time_end = time.perf_counter()
    print(f"MD simulation finished in {time_end - time_start:.2f} seconds")
    
    # 分析反应路径
    print("\n🔬 Analyzing reaction path...")
    path_info = detector.analyze_reaction_path()
    if path_info:
        print(f"   Total intermediates found: {path_info['total_intermediates']}")
        print(f"   Reaction steps: {path_info['reaction_steps']}")
        print(f"   Energy profile: {[f'{e:.3f}' for e in path_info['energy_profile']]}")
        
        # 保存反应路径分析
        with open("reaction_analysis/reaction_path.json", 'w') as f:
            json.dump(path_info, f, indent=2, default=str)
        
        print("   Reaction path analysis saved to reaction_analysis/reaction_path.json")
    else:
        print("   No reaction intermediates detected")


if __name__ == "__main__":
    main()



# 示例优化函数（已被后面的吸附优化函数替代）
# def objective(trial):
#     x = trial.suggest_float("x", 0, 1)
#     return x ** 2


def get_opt_energy(atoms, calculator, fmax=0.001, opt_mode: str = "normal"):    
    atoms.set_calculator(calculator)
    if opt_mode == "scale":
        opt1 = LBFGS(StrainFilter(atoms, mask=[1, 1, 1, 0, 0, 0]))
    elif opt_mode == "all":
        opt1 = LBFGS(ExpCellFilter(atoms))
    else:
        opt1 = LBFGS(atoms)
    opt1.run(fmax=fmax)
    return atoms.get_total_energy()



def already_slab(calculator):
    slab = read("surface.cif")
    slab.calc = calculator
    E_slab = get_opt_energy(slab, calculator, fmax=1e-4, opt_mode="normal")
    return slab, E_slab 
# 这些函数调用需要在有计算器的情况下使用
# slab, E_slab = already_slab(calculator)
# view_ngl(slab, representations=["ball+stick"])


def already_mol(filename, calculator):
    mol = read(filename)
    mol.calc = calculator  
    E_mol = get_opt_energy(mol, calculator, fmax=1e-4)
    return mol, E_mol


def get_all_intermediates(calculator):
    molecules = []
    for fname in [ "4-ketone.cif", "5-ketone.cif",]:
        mol, E_mol = already_mol(fname, calculator)
        molecules.append((mol, E_mol, fname))
    return molecules




def atoms_to_json(atoms):
    f = io.StringIO()
    write(f, atoms, format="json")
    return f.getvalue()


def json_to_atoms(atoms_str):
    return read(io.StringIO(atoms_str), format="json")

# 这些函数调用需要计算器参数，在实际使用时需要提供
# for mol, E_mol, name in get_all_intermediates(calculator):
#     print(f"🔍 Molecule: {name}")
#     mol_json_str = atoms_to_json(mol)
#     mol2 = json_to_atoms(mol_json_str)
#     
#     print(f"{mol_json_str=}")
#     view_ngl(mol2, representations=["ball+stick"]) 



def objective(trial):
    slab = json_to_atoms(trial.study.user_attrs["slab"])
    E_slab = trial.study.user_attrs["E_slab"]
    mol = json_to_atoms(trial.study.user_attrs["mol"])
    E_mol = trial.study.user_attrs["E_mol"]

    phi = 180. * trial.suggest_float("phi", -1, 1)
    theta = np.arccos(trial.suggest_float("theta", -1, 1)) * 180. / np.pi
    psi = 180. * trial.suggest_float("psi", -1, 1)
    x_pos = trial.suggest_float("x_pos", 0, 0.5)
    y_pos = trial.suggest_float("y_pos", 0, 0.5)
    z_hig = trial.suggest_float("z_hig", 2.0, 6.0)

    mol.euler_rotate(phi=phi, theta=theta, psi=psi)

    xy_position = np.matmul([x_pos, y_pos, 0], slab.cell)[:2]
    mol.translate([*xy_position, 0.0])

    max_slab_z = np.max(slab.get_positions()[:, 2])
    min_mol_z = np.min(mol.get_positions()[:, 2])
    shift_z = (max_slab_z + z_hig) - min_mol_z
    mol.translate([0, 0, shift_z])

    combined = slab + mol
    combined.calc = calculator

    E_slab_mol = get_opt_energy(combined, calculator, fmax=1e-3)
    trial.set_user_attr("structure", atoms_to_json(combined))

    return E_slab_mol - E_slab - E_mol



# 这些函数调用需要计算器参数，在实际使用时需要提供
# for mol, E_mol, name in get_all_intermediates(calculator):
#     print(f"\n🔍 Starting optimization for {name}...\n")
#     slab, E_slab = already_slab(calculator)
#     view(mol, viewer="ngl")
# 
#     study = optuna.create_study()
#     study.set_user_attr("slab", atoms_to_json(slab))
#     study.set_user_attr("E_slab", E_slab)
#     study.set_user_attr("mol", atoms_to_json(mol))
#     study.set_user_attr("E_mol", E_mol)
# 
#     study.optimize(objective, n_trials=300)
# 
#     print(f"    Best trial for {name} is #{study.best_trial.number}")
#     print(f"    Adsorption energy: {study.best_value:.6f} eV")
#     print("    Adsorption position:")
#     for key in ["phi", "theta", "psi", "x_pos", "y_pos", "z_hig"]:
#         print(f"        {key}: {study.best_params[key]}")
# 
#     output_dir = os.path.join("output", name)
#     os.makedirs(output_dir, exist_ok=True)
# 
#     # Save optimization plots
#     optuna.visualization.plot_optimization_history(study).write_html(
#         os.path.join(output_dir, "optimization_history.html"))
#     optuna.visualization.plot_slice(study).write_html(
#         os.path.join(output_dir, "optimization_slice.html"))
# 
#     # Save best trial structure
#     best_slab = json_to_atoms(study.best_trial.user_attrs["structure"])
#     view_ngl(best_slab, representations=["ball+stick"])
#     write(os.path.join(output_dir, "best_trial.cif"), best_slab)
# 
#     # Save all trials' structures and images
#     trial_data = []
#     n_trials = len(study.trials)
#     fig_rows = (n_trials // 10) + 1
#     fig, axes = plt.subplots(fig_rows, 10, figsize=(20, 2 * fig_rows))
# 
#     if fig_rows == 1:
#         axes = [axes]  # flatten if only one row
# 
#     for trial in study.trials:
#         slab = json_to_atoms(trial.user_attrs["structure"])
# 
#         # Save structure file
#         trial_cif = os.path.join(output_dir, f"{trial.number}.cif")
#         write(trial_cif, slab)
# 
#         # Save image
#         img_path = os.path.join(output_dir, f"{trial.number}.png")
#         write(img_path, slab, rotation="0x,0y,90z")
#         ax = axes[trial.number // 10][trial.number % 10]
#         ax.imshow(mpimg.imread(img_path))
#         ax.set_axis_off()
#         ax.set_title(trial.number)
# 
#         # Collect trial data
#         trial_data.append({
#             "trial": trial.number,
#             "energy": trial.value,
#             **trial.params
#         })
# 
#     # Save energy data CSV
#     df = pd.DataFrame(trial_data)
#     df.to_csv(os.path.join(output_dir, "data.csv"), index=False)
# 
#     # Save grid plot of all trials
#     plt.tight_layout()
#     plt.savefig(os.path.join(output_dir, "trial_grid.png"))
#     plt.show()
# 
#     # Show all trial structures in viewer
#     slabs = [json_to_atoms(trial.user_attrs["structure"]) for trial in study.trials]
#     view_ngl(slabs, representations=["ball+stick"], replace_structure=True)


def analyze_oxidation_reaction(atoms, calculator, temperature=300, timestep=0.5, steps=5000):
    """专门分析氧化反应的主函数"""
    print("🔥 Starting oxidation reaction analysis...")
    
    # 设置计算器
    atoms.calc = calculator
    
    # 初始优化
    print("   Optimizing initial structure...")
    opt = LBFGS(atoms, logfile=None)
    opt.run(fmax=0.01)
    initial_energy = atoms.get_potential_energy()
    print(f"   Initial energy: {initial_energy:.4f} eV")
    
    # 分析初始键结构
    detector = ReactionDetector(atoms, calculator, bond_threshold=2.0)
    initial_bonds = detector.analyze_bonds(atoms)
    print(f"   Initial bonds: {len(initial_bonds)}")
    
    # 设置MD参数
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature, force_temp=True)
    Stationary(atoms)
    
    # 创建NVT动力学
    dyn = NPT(
        atoms,
        timestep * units.fs,
        temperature_K=temperature,
        externalstress=0,
        ttime=100 * units.fs,
        pfactor=None,
        trajectory="oxidation_trajectory.xyz",
        loginterval=100,
    )
    
    # 创建监控器
    detector, tracker = create_reaction_monitor(atoms, calculator, "oxidation_analysis")
    
    # 反应监控回调
    def oxidation_monitor():
        step = dyn.get_number_of_steps()
        is_intermediate, intermediate = monitor_reaction_step(atoms, step, detector, tracker)
        
        if is_intermediate:
            print(f"🔥 Oxidation intermediate at step {step}")
            for change in intermediate['bond_changes']:
                if 'O' in change['bond_type']:
                    print(f"   Oxygen-related change: {change['type']} {change['bond_type']}")
    
    # 附加监控
    dyn.attach(oxidation_monitor, interval=100)
    
    # 运行MD
    print("   Running MD simulation...")
    dyn.run(steps)
    
    # 分析结果
    print("\n🔥 Oxidation reaction analysis complete!")
    path_info = detector.analyze_reaction_path()
    
    if path_info:
        print(f"   Intermediates found: {path_info['total_intermediates']}")
        
        # 专门分析氧相关的键变化
        oxygen_changes = {}
        for bond_type, changes in path_info['bond_changes_summary'].items():
            if 'O' in bond_type:
                oxygen_changes[bond_type] = changes
        
        if oxygen_changes:
            print("   Oxygen-related bond changes:")
            for bond_type, changes in oxygen_changes.items():
                print(f"     {bond_type}: {changes}")
        
        # 保存详细分析
        oxidation_analysis = {
            'initial_energy': initial_energy,
            'initial_bonds': initial_bonds,
            'reaction_path': path_info,
            'oxygen_changes': oxygen_changes,
            'temperature': temperature,
            'steps': steps
        }
        
        with open("oxidation_analysis/oxidation_summary.json", 'w') as f:
            json.dump(oxidation_analysis, f, indent=2, default=str)
        
        print("   Detailed analysis saved to oxidation_analysis/oxidation_summary.json")
    
    return detector, tracker


def find_transition_states(intermediates, calculator):
    """寻找过渡态"""
    print("🎯 Searching for transition states...")
    
    if len(intermediates) < 2:
        print("   Need at least 2 intermediates to find transition states")
        return []
    
    transition_states = []
    
    for i in range(len(intermediates) - 1):
        initial = intermediates[i]['atoms']
        final = intermediates[i + 1]['atoms']
        
        print(f"   Analyzing transition between intermediate {i} and {i+1}")
        
        # 使用NEB寻找过渡态
        try:
            # 创建NEB路径
            images = [initial.copy() for _ in range(5)]
            images.append(final.copy())
            
            for image in images:
                image.calc = calculator
            
            # 使用NEB
            neb = NEB(images)
            neb.interpolate()
            
            # 优化NEB
            optimizer = BFGS(neb, logfile=None)
            optimizer.run(fmax=0.05)
            
            # 寻找能量最高点（过渡态）
            energies = [image.get_potential_energy() for image in images]
            max_energy_idx = np.argmax(energies)
            
            transition_state = {
                'atoms': images[max_energy_idx].copy(),
                'energy': energies[max_energy_idx],
                'path_energies': energies,
                'between_intermediates': (i, i+1)
            }
            
            transition_states.append(transition_state)
            
            # 保存过渡态
            write(f"transition_states/ts_{i}_{i+1}.cif", images[max_energy_idx])
            
            print(f"     Transition state found with energy {energies[max_energy_idx]:.4f} eV")
            
        except Exception as e:
            print(f"     Failed to find transition state: {e}")
    
    return transition_states


def visualize_reaction_network(intermediates, transition_states, output_file="reaction_network.html"):
    """可视化反应网络"""
    print("📊 Creating reaction network visualization...")
    
    import plotly.graph_objects as go
    import plotly.offline as pyo
    
    # 准备数据
    node_x = []
    node_y = []
    node_text = []
    node_energies = []
    
    # 添加中间体
    for i, intermediate in enumerate(intermediates):
        node_x.append(i * 2)
        node_y.append(0)
        node_text.append(f"Intermediate {i}<br>E: {intermediate['energy']:.3f} eV")
        node_energies.append(intermediate['energy'])
    
    # 添加过渡态
    for ts in transition_states:
        idx1, idx2 = ts['between_intermediates']
        ts_x = (idx1 + idx2) * 2 / 2
        ts_y = 1
        node_x.append(ts_x)
        node_y.append(ts_y)
        node_text.append(f"TS<br>E: {ts['energy']:.3f} eV")
        node_energies.append(ts['energy'])
    
    # 创建节点
    fig = go.Figure()
    
    # 中间体节点
    fig.add_trace(go.Scatter(
        x=node_x[:len(intermediates)],
        y=node_y[:len(intermediates)],
        mode='markers+text',
        marker=dict(size=20, color='blue'),
        text=[f"I{i}" for i in range(len(intermediates))],
        textposition="middle center",
        name="Intermediates"
    ))
    
    # 过渡态节点
    if transition_states:
        fig.add_trace(go.Scatter(
            x=node_x[len(intermediates):],
            y=node_y[len(intermediates):],
            mode='markers+text',
            marker=dict(size=15, color='red'),
            text=["TS"] * len(transition_states),
            textposition="middle center",
            name="Transition States"
        ))
    
    # 添加连接线
    for ts in transition_states:
        idx1, idx2 = ts['between_intermediates']
        fig.add_trace(go.Scatter(
            x=[idx1 * 2, (idx1 + idx2) * 2 / 2, idx2 * 2],
            y=[0, 1, 0],
            mode='lines',
            line=dict(color='gray', width=2),
            showlegend=False
        ))
    
    fig.update_layout(
        title="Oxidation Reaction Network",
        xaxis_title="Reaction Coordinate",
        yaxis_title="Energy Profile",
        showlegend=True,
        hovermode='closest'
    )
    
    # 保存HTML文件
    pyo.plot(fig, filename=output_file, auto_open=False)
    print(f"   Reaction network saved to {output_file}")
    
    return fig


def run_oxidation_analysis_example():
    """运行氧化反应分析示例"""
    print("🚀 Starting oxidation reaction analysis example...")
    
    # 假设您有一个包含大分子和表面的结构
    # 这里使用示例结构，您需要替换为您的实际结构文件
    
    try:
        # 读取结构（请替换为您的实际结构文件）
        atoms = read("surface_with_molecule.cif")  # 请替换为您的结构文件
        
        # 设置计算器
        estimator = Estimator(model_version="v8.0.0", calc_mode=EstimatorCalcMode.CRYSTAL_U0_PLUS_D3)
        calculator = ASECalculator(estimator)
        
        # 运行氧化反应分析
        detector, tracker = analyze_oxidation_reaction(
            atoms, 
            calculator, 
            temperature=400,  # 较高温度促进反应
            timestep=0.5, 
            steps=10000
        )
        
        # 寻找过渡态
        if detector.intermediates:
            os.makedirs("transition_states", exist_ok=True)
            transition_states = find_transition_states(detector.intermediates, calculator)
            
            # 创建反应网络可视化
            if transition_states:
                fig = visualize_reaction_network(detector.intermediates, transition_states)
                
                print("\n✅ Oxidation analysis complete!")
                print(f"   Found {len(detector.intermediates)} intermediates")
                print(f"   Found {len(transition_states)} transition states")
                print("   Results saved in oxidation_analysis/ and transition_states/ directories")
            else:
                print("\n⚠️  No transition states found")
        else:
            print("\n⚠️  No reaction intermediates detected")
            
    except FileNotFoundError:
        print("❌ Structure file not found. Please provide your structure file.")
        print("   Expected file: surface_with_molecule.cif")
    except Exception as e:
        print(f"❌ Error during analysis: {e}")


def run_carbon_optimization_analysis(structure_file, temperature=300, md_steps=100000, opt_trials=100):
    """
    运行C原子优化分析的主函数
    
    Args:
        structure_file: 反应中间体结构文件路径
        temperature: MD模拟温度 (K)
        md_steps: MD模拟步数
        opt_trials: Optuna优化试验次数
    """
    print("🔬 Starting Carbon Optimization Analysis")
    print("=" * 50)
    
    try:
        # 1. 读取结构文件
        print(f"📖 Reading structure file: {structure_file}")
        atoms = read(structure_file)
        print(f"   Structure loaded: {len(atoms)} atoms")
        
        # 2. 设置计算器
        print("⚙️  Setting up calculator...")
        estimator = Estimator(model_version="v8.0.0", calc_mode=EstimatorCalcMode.CRYSTAL_U0_PLUS_D3)
        calculator = ASECalculator(estimator)
        atoms.calc = calculator
        
        # 3. 初始结构优化
        print("🔧 Initial structure optimization...")
        opt = LBFGS(atoms)
        opt.run(fmax=0.01)
        initial_energy = atoms.get_potential_energy()
        print(f"   Initial energy: {initial_energy:.4f} eV")
        
        # 4. 识别官能团和周围的C原子
        print("\n🔍 Identifying functional groups and nearby carbons...")
        analyzer = FunctionalGroupAnalyzer(atoms, calculator)
        functional_groups = analyzer.identify_functional_groups()
        print(f"   Found {len(functional_groups)} functional groups:")
        for group in functional_groups:
            print(f"     - {group['name']}: atoms {group['atoms']}")
        
        carbon_atoms = analyzer.find_carbons_near_functional_groups(max_distance=3.0)
        print(f"   Found {len(carbon_atoms)} carbon atoms near functional groups:")
        for i, carbon in enumerate(carbon_atoms[:5]):  # 只显示前5个
            print(f"     - Carbon {carbon['index']}: distance {carbon['distance_to_group']:.2f} Å to {carbon['functional_group']}")
        
        if len(carbon_atoms) == 0:
            print("⚠️  No carbon atoms found near functional groups!")
            return None
        
        # 5. MD模拟寻找最稳定结构
        print(f"\n🧪 Running MD simulation to find stable structure...")
        md_searcher = MDStabilitySearcher(
            atoms, 
            calculator, 
            temperature=temperature, 
            steps=md_steps, 
            sample_interval=100
        )
        
        most_stable = md_searcher.find_most_stable_structure()
        print(f"   Most stable structure energy: {most_stable['energy']:.4f} eV")
        
        # 使用最稳定的结构进行后续优化
        atoms = most_stable['atoms'].copy()
        
        # 6. 从近到远逐步优化C原子
        print(f"\n🎯 Starting progressive carbon optimization...")
        optimizer = ProgressiveCarbonOptimizer(
            atoms, 
            calculator, 
            carbon_atoms, 
            functional_groups
        )
        
        optimization_history = optimizer.optimize_carbons_progressively(n_trials=opt_trials)
        
        # 7. 保存结果
        print(f"\n💾 Saving results...")
        output_dir = "carbon_optimization_results"
        os.makedirs(output_dir, exist_ok=True)
        
        # 保存最终优化结构
        final_energy = atoms.get_potential_energy()
        write(f"{output_dir}/final_optimized_structure.cif", atoms)
        
        # 保存优化历史
        optimization_summary = {
            'initial_energy': initial_energy,
            'md_stable_energy': most_stable['energy'],
            'final_energy': final_energy,
            'energy_improvement': initial_energy - final_energy,
            'functional_groups': functional_groups,
            'carbon_atoms': carbon_atoms,
            'optimization_history': optimization_history
        }
        
        with open(f"{output_dir}/optimization_summary.json", 'w') as f:
            json.dump(optimization_summary, f, indent=2, default=str)
        
        # 8. 生成分析报告
        print(f"\n📊 Analysis Summary:")
        print(f"   Initial energy: {initial_energy:.4f} eV")
        print(f"   MD stable energy: {most_stable['energy']:.4f} eV")
        print(f"   Final optimized energy: {final_energy:.4f} eV")
        print(f"   Total energy improvement: {initial_energy - final_energy:.4f} eV")
        print(f"   Carbon atoms optimized: {len(carbon_atoms)}")
        print(f"   Optimization trials per carbon: {opt_trials}")
        
        print(f"\n✅ Analysis complete!")
        print(f"   Results saved in: {output_dir}/")
        print(f"   Final structure: {output_dir}/final_optimized_structure.cif")
        
        return {
            'atoms': atoms,
            'optimization_history': optimization_history,
            'summary': optimization_summary
        }
        
    except FileNotFoundError:
        print(f"❌ Structure file not found: {structure_file}")
        return None
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == "__main__":
    # 可以选择运行不同的分析模式
    import sys
    
    if len(sys.argv) > 1:
        if sys.argv[1] == "--oxidation":
            # 运行氧化反应分析
            run_oxidation_analysis_example()
        elif sys.argv[1] == "--carbon-optimization":
            # 运行C原子优化分析
            if len(sys.argv) > 2:
                structure_file = sys.argv[2]
                # 可选参数
                temperature = int(sys.argv[3]) if len(sys.argv) > 3 else 300
                md_steps = int(sys.argv[4]) if len(sys.argv) > 4 else 100000
                opt_trials = int(sys.argv[5]) if len(sys.argv) > 5 else 100
                
                run_carbon_optimization_analysis(
                    structure_file=structure_file,
                    temperature=temperature,
                    md_steps=md_steps,
                    opt_trials=opt_trials
                )
            else:
                print("❌ Please provide structure file path")
                print("Usage: python search.py --carbon-optimization <structure_file> [temperature] [md_steps] [opt_trials]")
                print("Example: python search.py --carbon-optimization 2-nonanone.cif 400 100000 100")
        else:
            print("❌ Unknown option. Available options:")
            print("  --oxidation: Run oxidation reaction analysis")
            print("  --carbon-optimization <file>: Run carbon optimization analysis")
            print("  (no option): Run standard MD simulation")
    else:
        # 运行标准MD程序
        main()