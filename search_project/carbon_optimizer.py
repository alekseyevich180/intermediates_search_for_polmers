#!/usr/bin/env python3
"""
C原子优化独立脚本

功能: 从近到远逐步优化C原子位置
输入: 结构文件 (支持表面+分子或纯分子)
输出: 优化后的结构、优化历史、能量变化
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
    from ase.constraints import FixAtoms
    from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
    from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
    import numpy as np
    import optuna
except ImportError as e:
    print(f"❌ 导入错误: {e}")
    print("请确保已安装所需的依赖包")
    sys.exit(1)


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
        study.set_user_attr("atoms", self._atoms_to_json(self.atoms))
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
    
    def _atoms_to_json(self, atoms):
        """将原子结构转换为JSON格式"""
        return {
            'positions': atoms.get_positions().tolist(),
            'symbols': atoms.get_chemical_symbols(),
            'cell': atoms.get_cell().tolist() if atoms.cell else None
        }


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


def identify_functional_groups(atoms):
    """识别分子中的官能团"""
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    
    functional_groups = []
    
    # 查找酮基 (C=O)
    for i, symbol in enumerate(symbols):
        if symbol == 'C':
            for j, other_symbol in enumerate(symbols):
                if other_symbol == 'O' and i != j:
                    distance = np.linalg.norm(positions[i] - positions[j])
                    if distance < 1.5:  # C=O键长
                        functional_groups.append({
                            'name': 'ketone',
                            'atoms': [i, j],
                            'center': (positions[i] + positions[j]) / 2
                        })
                        break
    
    return functional_groups


def find_carbons_near_functional_groups(atoms, functional_groups, max_distance=3.0):
    """找到官能团周围的C原子"""
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    
    nearby_carbons = []
    
    for group in functional_groups:
        group_center = group['center']
        
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
    return nearby_carbons


def run_carbon_optimization(structure_file, temperature=300, md_steps=100000, opt_trials=100):
    """运行C原子优化"""
    print(f"🔧 Running carbon optimization for: {structure_file}")
    
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
        
        # 识别官能团
        print("\n🔍 Identifying functional groups...")
        functional_groups = identify_functional_groups(atoms)
        print(f"   Found {len(functional_groups)} functional groups:")
        for group in functional_groups:
            print(f"     - {group['name']}: atoms {group['atoms']}")
        
        if not functional_groups:
            print("   ⚠️  No functional groups found, using all C atoms")
            # 如果没有找到官能团，使用所有C原子
            carbon_atoms = []
            for i, symbol in enumerate(symbols):
                if symbol == 'C':
                    carbon_atoms.append({
                        'index': i,
                        'distance_to_group': 0.0,
                        'functional_group': 'all_carbons',
                        'position': atoms.get_positions()[i].copy()
                    })
        else:
            # 找到周围的C原子
            print(f"\n⚛️  Finding carbon atoms near functional groups...")
            carbon_atoms = find_carbons_near_functional_groups(atoms, functional_groups)
            print(f"   Found {len(carbon_atoms)} carbon atoms:")
            for i, carbon in enumerate(carbon_atoms):
                print(f"     - Carbon {carbon['index']}: distance {carbon['distance_to_group']:.2f} Å")
        
        if not carbon_atoms:
            print("   ❌ No carbon atoms found for optimization")
            return None
        
        # 创建优化器
        optimizer = ProgressiveCarbonOptimizer(atoms, calculator, carbon_atoms, functional_groups)
        
        # 运行渐进式优化
        print(f"\n🔧 Starting progressive optimization...")
        start_time = time.time()
        optimization_history = optimizer.optimize_carbons_progressively(n_trials=opt_trials)
        end_time = time.time()
        
        # 获取最终能量
        final_energy = atoms.get_potential_energy()
        energy_improvement = initial_energy - final_energy
        
        print(f"\n✅ Optimization complete!")
        print(f"   Initial energy: {initial_energy:.4f} eV")
        print(f"   Final energy: {final_energy:.4f} eV")
        print(f"   Energy improvement: {energy_improvement:.4f} eV")
        print(f"   Optimization time: {end_time - start_time:.2f} seconds")
        
        # 保存结果
        output_dir = "carbon_optimization_results"
        os.makedirs(output_dir, exist_ok=True)
        
        # 保存优化后的结构
        write(f"{output_dir}/final_optimized_structure.cif", atoms)
        
        # 保存优化历史
        optimization_report = {
            'structure_file': structure_file,
            'initial_energy': initial_energy,
            'final_energy': final_energy,
            'energy_improvement': energy_improvement,
            'optimization_time': end_time - start_time,
            'functional_groups': functional_groups,
            'carbon_atoms': carbon_atoms,
            'optimization_history': optimization_history,
            'optimization_parameters': {
                'opt_trials': opt_trials,
                'has_surface': has_metal
            }
        }
        
        import json
        with open(f"{output_dir}/optimization_summary.json", 'w') as f:
            json.dump(optimization_report, f, indent=2, default=str)
        
        print(f"\n📁 Results saved in: {output_dir}/")
        print(f"   - final_optimized_structure.cif")
        print(f"   - optimization_summary.json")
        
        return {
            'final_atoms': atoms,
            'optimization_report': optimization_report,
            'energy_improvement': energy_improvement
        }
        
    except Exception as e:
        print(f"❌ Error during carbon optimization: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="Carbon Optimization Tool")
    parser.add_argument("structure_file", help="Structure file path (CIF, XYZ, etc.)")
    parser.add_argument("--opt_trials", type=int, default=100,
                       help="Number of optimization trials per carbon (default: 100)")
    parser.add_argument("--output_dir", type=str, default="carbon_optimization_results",
                       help="Output directory name (default: carbon_optimization_results)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.structure_file):
        print(f"❌ Structure file not found: {args.structure_file}")
        return
    
    print("🔧 Carbon Optimizer")
    print("=" * 50)
    print(f"📖 Structure file: {args.structure_file}")
    print(f"🔄 Optimization trials: {args.opt_trials}")
    print(f"📁 Output directory: {args.output_dir}")
    print()
    
    # 运行C原子优化
    result = run_carbon_optimization(
        args.structure_file,
        opt_trials=args.opt_trials
    )
    
    if result:
        print("\n🎉 Carbon optimization completed successfully!")
    else:
        print("\n❌ Carbon optimization failed!")


if __name__ == "__main__":
    main()
