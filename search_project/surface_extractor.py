#!/usr/bin/env python3
"""
表面分子提取独立脚本

功能: 从包含表面的结构中提取分子部分
输入: 表面+分子结构文件
输出: 提取的分子结构、分析报告、可视化信息
"""

import argparse
import os
import sys
from pathlib import Path

try:
    from ase.io import read, write
    from ase import Atoms
    import numpy as np
    import json
except ImportError as e:
    print(f"❌ 导入错误: {e}")
    print("请确保已安装所需的依赖包")
    sys.exit(1)


class SurfaceMoleculeExtractor:
    """从表面结构中提取分子部分"""
    
    def __init__(self, structure_file):
        self.structure_file = structure_file
        self.full_atoms = None
        self.surface_atoms = None
        self.molecule_atoms = None
        self.analysis_report = {}
        
    def load_structure(self):
        """加载结构文件"""
        print(f"📖 Loading structure: {self.structure_file}")
        
        try:
            self.full_atoms = read(self.structure_file)
            print(f"   ✅ Structure loaded: {len(self.full_atoms)} atoms")
            return True
        except Exception as e:
            print(f"   ❌ Error loading structure: {e}")
            return False
    
    def analyze_elements(self):
        """分析元素组成"""
        if self.full_atoms is None:
            return None
        
        symbols = self.full_atoms.get_chemical_symbols()
        unique_symbols = list(set(symbols))
        
        element_counts = {}
        for symbol in unique_symbols:
            element_counts[symbol] = symbols.count(symbol)
        
        print(f"   🧪 Element types: {unique_symbols}")
        print(f"   📊 Element counts: {element_counts}")
        
        return {
            'unique_symbols': unique_symbols,
            'element_counts': element_counts,
            'total_atoms': len(symbols)
        }
    
    def detect_metal_surface(self):
        """检测金属表面"""
        if self.full_atoms is None:
            return False
        
        symbols = self.full_atoms.get_chemical_symbols()
        metal_elements = ['Pt', 'Pd', 'Au', 'Ag', 'Cu', 'Ni', 'Fe', 'Ti', 'Al', 'Zn', 'Ir', 'Ru']
        
        metal_atoms = []
        for i, symbol in enumerate(symbols):
            if symbol in metal_elements:
                metal_atoms.append(i)
        
        has_metal = len(metal_atoms) > 0
        
        if has_metal:
            print(f"   🏗️  Metal surface detected: {len(metal_atoms)} metal atoms")
            for metal in metal_elements:
                count = symbols.count(metal)
                if count > 0:
                    print(f"     - {metal}: {count} atoms")
        else:
            print("   📝 No metal surface detected")
        
        return has_metal, metal_atoms
    
    def analyze_z_distribution(self):
        """分析Z坐标分布"""
        if self.full_atoms is None:
            return None
        
        positions = self.full_atoms.get_positions()
        z_coords = positions[:, 2]
        
        z_min = float(z_coords.min())
        z_max = float(z_coords.max())
        z_range = z_max - z_min
        z_mean = float(z_coords.mean())
        z_std = float(z_coords.std())
        
        print(f"   📏 Z coordinate analysis:")
        print(f"     - Range: {z_min:.2f} to {z_max:.2f} Å")
        print(f"     - Mean: {z_mean:.2f} Å")
        print(f"     - Std: {z_std:.2f} Å")
        
        return {
            'z_min': z_min,
            'z_max': z_max,
            'z_range': z_range,
            'z_mean': z_mean,
            'z_std': z_std,
            'z_coords': z_coords.tolist()
        }
    
    def identify_layers(self, z_analysis, surface_threshold=0.3, molecule_threshold=0.4):
        """识别表面层和分子层"""
        if z_analysis is None:
            return None
        
        z_coords = np.array(z_analysis['z_coords'])
        z_min = z_analysis['z_min']
        z_range = z_analysis['z_range']
        
        # 计算阈值
        surface_z = z_min + z_range * surface_threshold
        molecule_z = z_min + z_range * molecule_threshold
        
        # 识别各层原子
        surface_indices = np.where(z_coords < surface_z)[0]
        middle_indices = np.where((z_coords >= surface_z) & (z_coords <= molecule_z))[0]
        molecule_indices = np.where(z_coords > molecule_z)[0]
        
        print(f"   🔍 Layer identification:")
        print(f"     - Surface layer (z < {surface_z:.2f}): {len(surface_indices)} atoms")
        print(f"     - Middle layer ({surface_z:.2f} ≤ z ≤ {molecule_z:.2f}): {len(middle_indices)} atoms")
        print(f"     - Molecule layer (z > {molecule_z:.2f}): {len(molecule_indices)} atoms")
        
        return {
            'surface_indices': surface_indices.tolist(),
            'middle_indices': middle_indices.tolist(),
            'molecule_indices': molecule_indices.tolist(),
            'surface_threshold': surface_z,
            'molecule_threshold': molecule_z
        }
    
    def analyze_molecule_composition(self, layer_info):
        """分析分子层的组成"""
        if layer_info is None or len(layer_info['molecule_indices']) == 0:
            return None
        
        symbols = self.full_atoms.get_chemical_symbols()
        molecule_indices = layer_info['molecule_indices']
        
        molecule_symbols = [symbols[i] for i in molecule_indices]
        molecule_elements = list(set(molecule_symbols))
        
        molecule_counts = {}
        for element in molecule_elements:
            molecule_counts[element] = molecule_symbols.count(element)
        
        # 检查是否包含有机元素
        organic_elements = ['C', 'H', 'O', 'N', 'S', 'P']
        has_organic = any(elem in molecule_elements for elem in organic_elements)
        
        print(f"   🧪 Molecule layer analysis:")
        print(f"     - Elements: {molecule_elements}")
        print(f"     - Counts: {molecule_counts}")
        print(f"     - Has organic elements: {has_organic}")
        
        return {
            'elements': molecule_elements,
            'element_counts': molecule_counts,
            'has_organic': has_organic,
            'organic_elements': [elem for elem in molecule_elements if elem in organic_elements],
            'atom_indices': molecule_indices
        }
    
    def extract_molecule(self, layer_info, molecule_composition):
        """提取分子部分"""
        if (layer_info is None or molecule_composition is None or 
            not molecule_composition['has_organic']):
            print("   ⚠️  Cannot extract molecule: no organic elements found")
            return None
        
        molecule_indices = layer_info['molecule_indices']
        symbols = self.full_atoms.get_chemical_symbols()
        positions = self.full_atoms.get_positions()
        
        # 提取分子原子
        molecule_symbols = [symbols[i] for i in molecule_indices]
        molecule_positions = positions[molecule_indices]
        
        # 创建分子结构
        self.molecule_atoms = Atoms(symbols=molecule_symbols, positions=molecule_positions)
        
        print(f"   ✅ Molecule extracted: {len(self.molecule_atoms)} atoms")
        
        return self.molecule_atoms
    
    def extract_surface(self, layer_info):
        """提取表面部分"""
        if layer_info is None:
            return None
        
        surface_indices = layer_info['surface_indices']
        symbols = self.full_atoms.get_chemical_symbols()
        positions = self.full_atoms.get_positions()
        
        # 提取表面原子
        surface_symbols = [symbols[i] for i in surface_indices]
        surface_positions = positions[surface_indices]
        
        # 创建表面结构
        self.surface_atoms = Atoms(symbols=surface_symbols, positions=surface_positions)
        
        print(f"   ✅ Surface extracted: {len(self.surface_atoms)} atoms")
        
        return self.surface_atoms
    
    def save_extracted_structures(self, output_dir):
        """保存提取的结构"""
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        saved_files = []
        
        # 保存完整结构
        if self.full_atoms is not None:
            full_file = f"{output_dir}/full_structure.cif"
            write(full_file, self.full_atoms)
            saved_files.append(full_file)
            print(f"   💾 Full structure saved: {full_file}")
        
        # 保存表面结构
        if self.surface_atoms is not None:
            surface_file = f"{output_dir}/surface_only.cif"
            write(surface_file, self.surface_atoms)
            saved_files.append(surface_file)
            print(f"   💾 Surface saved: {surface_file}")
        
        # 保存分子结构
        if self.molecule_atoms is not None:
            molecule_file = f"{output_dir}/molecule_only.cif"
            write(molecule_file, self.molecule_atoms)
            saved_files.append(molecule_file)
            print(f"   💾 Molecule saved: {molecule_file}")
        
        return saved_files
    
    def generate_analysis_report(self):
        """生成分析报告"""
        self.analysis_report = {
            'structure_file': self.structure_file,
            'extraction_successful': self.molecule_atoms is not None,
            'total_atoms': len(self.full_atoms) if self.full_atoms else 0,
            'surface_atoms': len(self.surface_atoms) if self.surface_atoms else 0,
            'molecule_atoms': len(self.molecule_atoms) if self.molecule_atoms else 0,
            'timestamp': str(time.time()) if 'time' in globals() else 'unknown'
        }
        
        return self.analysis_report
    
    def run_extraction(self, surface_threshold=0.3, molecule_threshold=0.4):
        """运行完整的提取过程"""
        print("🔍 Surface Molecule Extractor")
        print("=" * 50)
        
        # 1. 加载结构
        if not self.load_structure():
            return False
        
        # 2. 分析元素组成
        element_analysis = self.analyze_elements()
        
        # 3. 检测金属表面
        has_metal, metal_atoms = self.detect_metal_surface()
        
        # 4. 分析Z坐标分布
        z_analysis = self.analyze_z_distribution()
        
        # 5. 识别各层
        layer_info = self.identify_layers(z_analysis, surface_threshold, molecule_threshold)
        
        # 6. 分析分子组成
        molecule_composition = self.analyze_molecule_composition(layer_info)
        
        # 7. 提取分子
        molecule_atoms = self.extract_molecule(layer_info, molecule_composition)
        
        # 8. 提取表面
        surface_atoms = self.extract_surface(layer_info)
        
        # 9. 生成报告
        report = self.generate_analysis_report()
        
        # 添加详细信息到报告
        report.update({
            'element_analysis': element_analysis,
            'has_metal_surface': has_metal,
            'metal_atoms_count': len(metal_atoms),
            'z_analysis': z_analysis,
            'layer_info': layer_info,
            'molecule_composition': molecule_composition,
            'extraction_parameters': {
                'surface_threshold': surface_threshold,
                'molecule_threshold': molecule_threshold
            }
        })
        
        self.analysis_report = report
        
        print(f"\n✅ Extraction complete!")
        print(f"   📊 Total atoms: {report['total_atoms']}")
        print(f"   🏗️  Surface atoms: {report['surface_atoms']}")
        print(f"   🧪 Molecule atoms: {report['molecule_atoms']}")
        print(f"   ✅ Extraction successful: {report['extraction_successful']}")
        
        return True


def run_surface_extraction(structure_file, output_dir="surface_extraction_results", 
                          surface_threshold=0.3, molecule_threshold=0.4):
    """运行表面分子提取"""
    print(f"🔍 Running surface extraction for: {structure_file}")
    
    try:
        # 创建提取器
        extractor = SurfaceMoleculeExtractor(structure_file)
        
        # 运行提取
        success = extractor.run_extraction(surface_threshold, molecule_threshold)
        
        if not success:
            print("❌ Extraction failed!")
            return None
        
        # 保存提取的结构
        saved_files = extractor.save_extracted_structures(output_dir)
        
        # 保存分析报告
        report_file = f"{output_dir}/extraction_analysis.json"
        with open(report_file, 'w') as f:
            json.dump(extractor.analysis_report, f, indent=2, default=str)
        
        print(f"\n📁 Results saved in: {output_dir}/")
        for file in saved_files:
            print(f"   - {os.path.basename(file)}")
        print(f"   - extraction_analysis.json")
        
        return {
            'extractor': extractor,
            'analysis_report': extractor.analysis_report,
            'saved_files': saved_files
        }
        
    except Exception as e:
        print(f"❌ Error during extraction: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="Surface Molecule Extractor Tool")
    parser.add_argument("structure_file", help="Surface+molecule structure file path")
    parser.add_argument("--output_dir", type=str, default="surface_extraction_results",
                       help="Output directory name (default: surface_extraction_results)")
    parser.add_argument("--surface_threshold", type=float, default=0.3,
                       help="Surface layer threshold (default: 0.3)")
    parser.add_argument("--molecule_threshold", type=float, default=0.4,
                       help="Molecule layer threshold (default: 0.4)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.structure_file):
        print(f"❌ Structure file not found: {args.structure_file}")
        return
    
    print("🔍 Surface Molecule Extractor")
    print("=" * 50)
    print(f"📖 Structure file: {args.structure_file}")
    print(f"📁 Output directory: {args.output_dir}")
    print(f"🔍 Surface threshold: {args.surface_threshold}")
    print(f"🔍 Molecule threshold: {args.molecule_threshold}")
    print()
    
    # 运行表面分子提取
    result = run_surface_extraction(
        args.structure_file,
        output_dir=args.output_dir,
        surface_threshold=args.surface_threshold,
        molecule_threshold=args.molecule_threshold
    )
    
    if result:
        print("\n🎉 Surface extraction completed successfully!")
    else:
        print("\n❌ Surface extraction failed!")


if __name__ == "__main__":
    main()
