#!/usr/bin/env python3
"""
反应检测独立脚本

功能: 检测和分析氧化反应中的中间状态
输入: 结构文件 (支持表面+分子或纯分子)
输出: 反应中间体、键变化分析、反应路径
"""

import argparse
import os
import sys
import time
from pathlib import Path

try:
    from ase.io import read, write
    from ase import Atoms
    from ase.optimize import LBFGS
    from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
    from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
    import numpy as np
    import json
except ImportError as e:
    print(f"❌ 导入错误: {e}")
    print("请确保已安装所需的依赖包")
    sys.exit(1)


class ReactionDetector:
    """检测和分析氧化反应中的中间状态"""
    
    def __init__(self, atoms, calculator):
        self.atoms = atoms
        self.calculator = calculator
        self.bond_network = {}
        self.reaction_intermediates = []
        self.bond_changes = []
        
    def analyze_bonds(self):
        """分析当前结构的化学键"""
        positions = self.atoms.get_positions()
        symbols = self.atoms.get_chemical_symbols()
        
        bonds = {}
        bond_distances = {}
        
        for i in range(len(self.atoms)):
            bonds[i] = []
            for j in range(i+1, len(self.atoms)):
                distance = np.linalg.norm(positions[i] - positions[j])
                if self._is_bond(i, j, distance, symbols):
                    bonds[i].append(j)
                    bonds[j].append(i)
                    bond_distances[(i, j)] = distance
                    bond_distances[(j, i)] = distance
        
        self.bond_network = bonds
        return bonds, bond_distances
    
    def _is_bond(self, i, j, distance, symbols):
        """判断两个原子是否成键"""
        bond_lengths = {
            ('C', 'C'): 1.8,
            ('C', 'O'): 1.7,
            ('C', 'N'): 1.7,
            ('C', 'H'): 1.3,
            ('O', 'H'): 1.2,
            ('N', 'H'): 1.2,
            ('O', 'O'): 1.6,
            ('N', 'N'): 1.6,
        }
        
        pair = tuple(sorted([symbols[i], symbols[j]]))
        max_bond_length = bond_lengths.get(pair, 2.5)
        
        return distance < max_bond_length
    
    def detect_bond_changes(self, previous_bonds, current_bonds):
        """检测键的变化"""
        changes = {
            'formed': [],
            'broken': [],
            'stretched': []
        }
        
        # 检测新形成的键
        for atom, bonds in current_bonds.items():
            for bond in bonds:
                if atom < bond:  # 避免重复
                    if (atom not in previous_bonds or 
                        bond not in previous_bonds[atom]):
                        changes['formed'].append((atom, bond))
        
        # 检测断裂的键
        for atom, bonds in previous_bonds.items():
            for bond in bonds:
                if atom < bond:  # 避免重复
                    if (atom not in current_bonds or 
                        bond not in current_bonds[atom]):
                        changes['broken'].append((atom, bond))
        
        return changes
    
    def is_reaction_intermediate(self, bond_changes, energy_threshold=0.1):
        """判断是否为反应中间体"""
        # 简化的判断标准
        has_bond_formation = len(bond_changes['formed']) > 0
        has_bond_breaking = len(bond_changes['broken']) > 0
        has_significant_changes = (len(bond_changes['formed']) + 
                                 len(bond_changes['broken'])) >= 2
        
        return has_bond_formation or has_bond_breaking or has_significant_changes
    
    def analyze_oxidation_reaction(self):
        """分析氧化反应过程"""
        print("🔍 Analyzing oxidation reaction...")
        
        # 初始结构分析
        initial_bonds, initial_distances = self.analyze_bonds()
        initial_energy = self.atoms.get_potential_energy()
        
        print(f"   Initial energy: {initial_energy:.4f} eV")
        print(f"   Initial bonds: {sum(len(bonds) for bonds in initial_bonds.values()) // 2}")
        
        # 存储初始状态
        self.reaction_intermediates.append({
            'step': 0,
            'atoms': self.atoms.copy(),
            'energy': initial_energy,
            'bonds': initial_bonds.copy(),
            'bond_distances': initial_distances.copy(),
            'type': 'initial'
        })
        
        # 模拟反应过程（通过多次优化和扰动）
        for step in range(1, 11):  # 模拟10个步骤
            print(f"   Step {step}: Analyzing reaction intermediate...")
            
            # 创建扰动后的结构
            perturbed_atoms = self.atoms.copy()
            self._apply_perturbation(perturbed_atoms)
            
            # 优化扰动后的结构
            perturbed_atoms.calc = self.calculator
            opt = LBFGS(perturbed_atoms)
            opt.run(fmax=0.01)
            
            # 分析新结构
            new_bonds, new_distances = self._analyze_bonds_for_atoms(perturbed_atoms)
            new_energy = perturbed_atoms.get_potential_energy()
            
            # 检测键变化
            bond_changes = self.detect_bond_changes(self.bond_network, new_bonds)
            
            # 判断是否为反应中间体
            if self.is_reaction_intermediate(bond_changes):
                print(f"     ✅ Reaction intermediate detected at step {step}")
                
                intermediate = {
                    'step': step,
                    'atoms': perturbed_atoms.copy(),
                    'energy': new_energy,
                    'bonds': new_bonds,
                    'bond_distances': new_distances,
                    'bond_changes': bond_changes,
                    'type': 'intermediate'
                }
                
                self.reaction_intermediates.append(intermediate)
                self.bond_changes.append(bond_changes)
                
                # 更新当前结构
                self.atoms = perturbed_atoms
                self.bond_network = new_bonds
            
            else:
                print(f"     ⚠️  No significant reaction detected at step {step}")
        
        return self.reaction_intermediates
    
    def _apply_perturbation(self, atoms):
        """对结构施加小的扰动"""
        positions = atoms.get_positions()
        
        # 添加小的随机位移
        noise = np.random.normal(0, 0.1, positions.shape)
        new_positions = positions + noise
        atoms.set_positions(new_positions)
    
    def _analyze_bonds_for_atoms(self, atoms):
        """分析给定原子结构的化学键"""
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        bonds = {}
        bond_distances = {}
        
        for i in range(len(atoms)):
            bonds[i] = []
            for j in range(i+1, len(atoms)):
                distance = np.linalg.norm(positions[i] - positions[j])
                if self._is_bond(i, j, distance, symbols):
                    bonds[i].append(j)
                    bonds[j].append(i)
                    bond_distances[(i, j)] = distance
                    bond_distances[(j, i)] = distance
        
        return bonds, bond_distances
    
    def generate_reaction_report(self):
        """生成反应分析报告"""
        if not self.reaction_intermediates:
            return None
        
        report = {
            'total_intermediates': len(self.reaction_intermediates),
            'initial_energy': self.reaction_intermediates[0]['energy'],
            'final_energy': self.reaction_intermediates[-1]['energy'],
            'energy_change': (self.reaction_intermediates[-1]['energy'] - 
                            self.reaction_intermediates[0]['energy']),
            'intermediates': []
        }
        
        for intermediate in self.reaction_intermediates:
            intermediate_info = {
                'step': intermediate['step'],
                'energy': intermediate['energy'],
                'bond_count': sum(len(bonds) for bonds in intermediate['bonds'].values()) // 2,
                'type': intermediate['type']
            }
            
            if 'bond_changes' in intermediate:
                intermediate_info['bond_changes'] = intermediate['bond_changes']
            
            report['intermediates'].append(intermediate_info)
        
        return report


def extract_molecule_from_surface(structure_file):
    """从包含表面的结构中提取分子部分"""
    print(f"🔍 Extracting molecule from surface structure: {structure_file}")
    
    try:
        # 读取完整结构
        full_atoms = read(structure_file)
        print(f"   Total atoms: {len(full_atoms)}")
        
        # 分析原子类型
        symbols = full_atoms.get_chemical_symbols()
        unique_symbols = list(set(symbols))
        
        # 识别表面和分子层
        positions = full_atoms.get_positions()
        z_coords = positions[:, 2]
        z_min, z_max = z_coords.min(), z_coords.max()
        z_range = z_max - z_min
        
        # 识别表面层 (底部30%)
        surface_threshold = z_min + z_range * 0.3
        surface_indices = np.where(z_coords < surface_threshold)[0]
        
        # 识别分子层 (40%以上)
        molecule_threshold = z_min + z_range * 0.4
        molecule_indices = np.where(z_coords > molecule_threshold)[0]
        
        print(f"   Surface atoms: {len(surface_indices)}")
        print(f"   Molecule atoms: {len(molecule_indices)}")
        
        # 检查分子部分是否包含有机元素
        if len(molecule_indices) > 0:
            molecule_symbols = [symbols[i] for i in molecule_indices]
            organic_elements = ['C', 'H', 'O', 'N', 'S', 'P']
            has_organic = any(elem in molecule_symbols for elem in organic_elements)
            
            if has_organic:
                print("   ✅ Organic molecule detected")
                
                # 提取分子部分
                molecule_positions = positions[molecule_indices]
                molecule_atoms = Atoms(symbols=molecule_symbols, positions=molecule_positions)
                
                # 保存分子部分
                molecule_file = structure_file.replace('.cif', '_molecule_only.cif')
                write(molecule_file, molecule_atoms)
                print(f"   💾 Molecule saved as: {molecule_file}")
                
                return molecule_file, molecule_atoms
            else:
                print("   ⚠️  No organic molecule detected")
                return structure_file, full_atoms
        else:
            print("   ⚠️  No molecule layer detected, using full structure")
            return structure_file, full_atoms
            
    except Exception as e:
        print(f"   ❌ Error extracting molecule: {e}")
        return structure_file, read(structure_file)


def run_reaction_detection(structure_file):
    """运行反应检测分析"""
    print(f"🔍 Running reaction detection for: {structure_file}")
    
    try:
        # 读取结构
        atoms = read(structure_file)
        print(f"   📊 Structure loaded: {len(atoms)} atoms")
        
        # 检查是否包含表面
        symbols = atoms.get_chemical_symbols()
        metal_elements = ['Pt', 'Pd', 'Au', 'Ag', 'Cu', 'Ni', 'Fe', 'Ti', 'Al', 'Zn', 'Ir', 'Ru']
        has_metal = any(metal in symbols for metal in metal_elements)
        
        if has_metal and len(atoms) > 50:
            print("   🔍 Metal surface detected, extracting molecule...")
            molecule_file, atoms = extract_molecule_from_surface(structure_file)
            print(f"   Using extracted molecule: {len(atoms)} atoms")
        else:
            print("   📝 Using structure as-is")
        
        # 设置计算器
        estimator = Estimator(model_version="v8.0.0", calc_mode=EstimatorCalcMode.CRYSTAL_U0_PLUS_D3)
        calculator = ASECalculator(estimator)
        atoms.calc = calculator
        
        # 初始结构优化
        print("\n🔧 Initial structure optimization...")
        opt = LBFGS(atoms)
        opt.run(fmax=0.01)
        initial_energy = atoms.get_potential_energy()
        print(f"   Initial energy: {initial_energy:.4f} eV")
        
        # 创建反应检测器
        detector = ReactionDetector(atoms, calculator)
        
        # 分析氧化反应
        print(f"\n🔍 Analyzing oxidation reaction...")
        start_time = time.time()
        intermediates = detector.analyze_oxidation_reaction()
        end_time = time.time()
        
        # 生成反应报告
        reaction_report = detector.generate_reaction_report()
        
        print(f"\n✅ Reaction detection complete!")
        print(f"   Analysis time: {end_time - start_time:.2f} seconds")
        print(f"   Total intermediates found: {len(intermediates)}")
        
        if reaction_report:
            print(f"   Initial energy: {reaction_report['initial_energy']:.4f} eV")
            print(f"   Final energy: {reaction_report['final_energy']:.4f} eV")
            print(f"   Energy change: {reaction_report['energy_change']:.4f} eV")
        
        # 保存结果
        output_dir = "reaction_detection_results"
        os.makedirs(output_dir, exist_ok=True)
        
        # 保存所有中间体
        for i, intermediate in enumerate(intermediates):
            filename = f"{output_dir}/intermediate_{intermediate['step']:02d}.cif"
            write(filename, intermediate['atoms'])
        
        # 保存反应报告
        if reaction_report:
            with open(f"{output_dir}/reaction_analysis.json", 'w') as f:
                json.dump(reaction_report, f, indent=2, default=str)
        
        # 保存键变化分析
        if detector.bond_changes:
            with open(f"{output_dir}/bond_changes.json", 'w') as f:
                json.dump(detector.bond_changes, f, indent=2, default=str)
        
        print(f"\n📁 Results saved in: {output_dir}/")
        print(f"   - intermediate_XX.cif (intermediate structures)")
        print(f"   - reaction_analysis.json (reaction report)")
        print(f"   - bond_changes.json (bond change analysis)")
        
        return {
            'intermediates': intermediates,
            'reaction_report': reaction_report,
            'bond_changes': detector.bond_changes
        }
        
    except Exception as e:
        print(f"❌ Error during reaction detection: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="Reaction Detection Tool")
    parser.add_argument("structure_file", help="Structure file path (CIF, XYZ, etc.)")
    parser.add_argument("--output_dir", type=str, default="reaction_detection_results",
                       help="Output directory name (default: reaction_detection_results)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.structure_file):
        print(f"❌ Structure file not found: {args.structure_file}")
        return
    
    print("🔍 Reaction Detector")
    print("=" * 50)
    print(f"📖 Structure file: {args.structure_file}")
    print(f"📁 Output directory: {args.output_dir}")
    print()
    
    # 运行反应检测
    result = run_reaction_detection(args.structure_file)
    
    if result:
        print("\n🎉 Reaction detection completed successfully!")
    else:
        print("\n❌ Reaction detection failed!")


if __name__ == "__main__":
    main()
