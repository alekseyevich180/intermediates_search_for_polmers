#!/usr/bin/env python3
"""
官能团分析独立脚本

功能: 识别分子中的官能团和周围的C原子
输入: 结构文件 (支持表面+分子或纯分子)
输出: 官能团信息、C原子位置、分析报告
"""

import argparse
import os
import sys
from pathlib import Path

try:
    from ase.io import read, write
    from ase import Atoms
    import numpy as np
    from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
    from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
except ImportError as e:
    print(f"❌ 导入错误: {e}")
    print("请确保已安装所需的依赖包")
    sys.exit(1)


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
        print(f"   Element types: {unique_symbols}")
        
        # 统计各元素数量
        element_counts = {}
        for symbol in unique_symbols:
            element_counts[symbol] = symbols.count(symbol)
        print(f"   Element counts: {element_counts}")
        
        # 识别表面和分子层
        positions = full_atoms.get_positions()
        z_coords = positions[:, 2]
        z_min, z_max = z_coords.min(), z_coords.max()
        z_range = z_max - z_min
        
        print(f"   Z coordinate range: {z_min:.2f} to {z_max:.2f} Å")
        
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


def analyze_functional_groups(structure_file, max_distance=3.0):
    """分析结构中的官能团和C原子"""
    print(f"🔍 Analyzing functional groups in: {structure_file}")
    
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
        
        # 创建分析器
        analyzer = FunctionalGroupAnalyzer(atoms, calculator)
        
        # 识别官能团
        print("\n🔍 Identifying functional groups...")
        functional_groups = analyzer.identify_functional_groups()
        print(f"   Found {len(functional_groups)} functional groups:")
        for group in functional_groups:
            print(f"     - {group['name']}: atoms {group['atoms']}")
            print(f"       Center: [{group['center'][0]:.2f}, {group['center'][1]:.2f}, {group['center'][2]:.2f}]")
        
        # 找到周围的C原子
        print(f"\n⚛️  Finding carbon atoms near functional groups (max distance: {max_distance} Å)...")
        carbon_atoms = analyzer.find_carbons_near_functional_groups(max_distance)
        print(f"   Found {len(carbon_atoms)} carbon atoms:")
        for i, carbon in enumerate(carbon_atoms):
            print(f"     - Carbon {carbon['index']}: distance {carbon['distance_to_group']:.2f} Å to {carbon['functional_group']}")
            print(f"       Position: [{carbon['position'][0]:.2f}, {carbon['position'][1]:.2f}, {carbon['position'][2]:.2f}]")
        
        # 构建键网络
        print(f"\n🔗 Building bond network...")
        bond_network = analyzer.build_bond_network()
        total_bonds = sum(len(bonds) for bonds in bond_network.values()) // 2
        print(f"   Total bonds detected: {total_bonds}")
        
        # 保存结果
        output_dir = "functional_group_analysis"
        os.makedirs(output_dir, exist_ok=True)
        
        # 保存分析报告
        analysis_report = {
            'structure_file': structure_file,
            'total_atoms': len(atoms),
            'functional_groups': functional_groups,
            'carbon_atoms': carbon_atoms,
            'bond_network_size': total_bonds,
            'analysis_parameters': {
                'max_distance': max_distance,
                'has_surface': has_metal
            }
        }
        
        import json
        with open(f"{output_dir}/functional_group_analysis.json", 'w') as f:
            json.dump(analysis_report, f, indent=2, default=str)
        
        # 保存分子结构（如果有提取的话）
        if has_metal:
            write(f"{output_dir}/extracted_molecule.cif", atoms)
            print(f"   💾 Extracted molecule saved as: {output_dir}/extracted_molecule.cif")
        
        print(f"\n✅ Analysis complete!")
        print(f"   📁 Results saved in: {output_dir}/")
        print(f"   📊 Functional groups: {len(functional_groups)}")
        print(f"   ⚛️  Carbon atoms: {len(carbon_atoms)}")
        print(f"   🔗 Total bonds: {total_bonds}")
        
        return analysis_report
        
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="Functional Group Analysis Tool")
    parser.add_argument("structure_file", help="Structure file path (CIF, XYZ, etc.)")
    parser.add_argument("--max_distance", type=float, default=3.0, 
                       help="Maximum distance to search for carbon atoms (default: 3.0 Å)")
    parser.add_argument("--output_dir", type=str, default="functional_group_analysis",
                       help="Output directory name (default: functional_group_analysis)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.structure_file):
        print(f"❌ Structure file not found: {args.structure_file}")
        return
    
    print("🔍 Functional Group Analyzer")
    print("=" * 50)
    print(f"📖 Structure file: {args.structure_file}")
    print(f"🔍 Max distance: {args.max_distance} Å")
    print(f"📁 Output directory: {args.output_dir}")
    print()
    
    # 运行分析
    result = analyze_functional_groups(args.structure_file, args.max_distance)
    
    if result:
        print("\n🎉 Analysis completed successfully!")
    else:
        print("\n❌ Analysis failed!")


if __name__ == "__main__":
    main()
