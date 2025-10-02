#!/usr/bin/env python3
"""
表面+分子结构分析脚本

专门用于处理包含表面的结构文件，自动识别和提取分子部分进行C原子优化
"""

import os
import sys
import subprocess
from pathlib import Path

def analyze_surface_molecule_structure(structure_file):
    """分析表面+分子结构"""
    print(f"🔍 Analyzing surface+molecule structure: {structure_file}")
    
    try:
        from ase.io import read
        from ase import Atoms
        import numpy as np
        
        # 读取结构
        full_atoms = read(structure_file)
        print(f"   📊 Total atoms: {len(full_atoms)}")
        
        # 分析元素组成
        symbols = full_atoms.get_chemical_symbols()
        unique_symbols = list(set(symbols))
        element_counts = {}
        for symbol in unique_symbols:
            element_counts[symbol] = symbols.count(symbol)
        
        print(f"   🧪 Element types: {unique_symbols}")
        print(f"   📈 Element counts: {element_counts}")
        
        # 按z坐标分层分析
        positions = full_atoms.get_positions()
        z_coords = positions[:, 2]
        z_min, z_max = z_coords.min(), z_coords.max()
        z_range = z_max - z_min
        
        print(f"   📏 Z coordinate range: {z_min:.2f} to {z_max:.2f} Å (range: {z_range:.2f} Å)")
        
        # 分层识别
        surface_threshold = z_min + z_range * 0.3
        molecule_threshold = z_min + z_range * 0.4
        
        surface_indices = np.where(z_coords < surface_threshold)[0]
        molecule_indices = np.where(z_coords > molecule_threshold)[0]
        middle_indices = np.where((z_coords >= surface_threshold) & (z_coords <= molecule_threshold))[0]
        
        print(f"   🏗️  Surface layer (z < {surface_threshold:.2f}): {len(surface_indices)} atoms")
        print(f"   🧪 Molecule layer (z > {molecule_threshold:.2f}): {len(molecule_indices)} atoms")
        print(f"   🔄 Middle layer: {len(middle_indices)} atoms")
        
        # 分析各层的元素组成
        if len(surface_indices) > 0:
            surface_symbols = [symbols[i] for i in surface_indices]
            surface_elements = list(set(surface_symbols))
            print(f"   🏗️  Surface elements: {surface_elements}")
        
        if len(molecule_indices) > 0:
            molecule_symbols = [symbols[i] for i in molecule_indices]
            molecule_elements = list(set(molecule_symbols))
            molecule_counts = {}
            for element in molecule_elements:
                molecule_counts[element] = molecule_symbols.count(element)
            
            print(f"   🧪 Molecule elements: {molecule_elements}")
            print(f"   📊 Molecule counts: {molecule_counts}")
            
            # 检查是否是有机分子
            organic_elements = ['C', 'H', 'O', 'N', 'S', 'P']
            has_organic = any(elem in molecule_elements for elem in organic_elements)
            
            if has_organic:
                print("   ✅ Organic molecule detected in molecule layer")
                
                # 提取分子部分
                molecule_positions = positions[molecule_indices]
                molecule_atoms = Atoms(symbols=molecule_symbols, positions=molecule_positions)
                
                # 保存分子部分
                molecule_file = structure_file.replace('.cif', '_extracted_molecule.cif')
                from ase.io import write
                write(molecule_file, molecule_atoms)
                
                print(f"   💾 Extracted molecule saved as: {molecule_file}")
                print(f"   📊 Extracted molecule atoms: {len(molecule_atoms)}")
                
                return {
                    'success': True,
                    'molecule_file': molecule_file,
                    'molecule_atoms': molecule_atoms,
                    'molecule_indices': molecule_indices,
                    'surface_atoms': len(surface_indices),
                    'molecule_atoms_count': len(molecule_atoms),
                    'total_atoms': len(full_atoms)
                }
            else:
                print("   ⚠️  No organic molecule detected in molecule layer")
                return {'success': False, 'reason': 'No organic molecule detected'}
        else:
            print("   ⚠️  No molecule layer detected")
            return {'success': False, 'reason': 'No molecule layer detected'}
            
    except Exception as e:
        print(f"   ❌ Error analyzing structure: {e}")
        return {'success': False, 'reason': str(e)}


def run_carbon_optimization_on_extracted_molecule(molecule_file, temperature=400, md_steps=10000, opt_trials=50):
    """对提取的分子运行C原子优化"""
    print(f"\n🎯 Running carbon optimization on extracted molecule: {molecule_file}")
    
    # 查找search.py文件
    current_dir = Path(__file__).parent
    search_py_path = None
    
    possible_paths = [
        current_dir / "search.py",
        current_dir.parent / "optuna+MD" / "search.py"
    ]
    
    for path in possible_paths:
        if path.exists():
            search_py_path = str(path.absolute())
            break
    
    if not search_py_path:
        print("❌ Cannot find search.py file!")
        return False
    
    # 构建命令
    cmd = [
        sys.executable,
        search_py_path,
        "--carbon-optimization",
        molecule_file,
        str(temperature),
        str(md_steps),
        str(opt_trials)
    ]
    
    print(f"🔧 Command: {' '.join(cmd)}")
    
    try:
        # 运行优化
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        
        if result.stdout:
            print("📊 Output:")
            print(result.stdout)
        
        if result.stderr:
            print("⚠️  Warnings/Errors:")
            print(result.stderr)
        
        if result.returncode == 0:
            print("✅ Carbon optimization completed successfully!")
            return True
        else:
            print(f"❌ Optimization failed with return code: {result.returncode}")
            return False
            
    except subprocess.TimeoutExpired:
        print("⏰ Optimization timed out after 1 hour")
        return False
    except Exception as e:
        print(f"❌ Error running optimization: {e}")
        return False


def create_sample_surface_molecule():
    """创建示例表面+分子结构"""
    print("🔧 Creating sample surface+molecule structure...")
    
    try:
        from ase import Atoms
        from ase.io import write
        import numpy as np
        
        # 创建金属表面 (Pt)
        surface_positions = []
        surface_symbols = []
        
        # 3x3 Pt表面，两层
        for i in range(3):
            for j in range(3):
                # 底层
                surface_positions.append([i*2.8, j*2.8, 0.0])
                surface_symbols.append('Pt')
                # 顶层
                surface_positions.append([i*2.8, j*2.8, 2.3])
                surface_symbols.append('Pt')
        
        # 创建有机分子 (2-nonanone) 在表面上
        molecule_positions = [
            # C链 (在表面上方的z位置)
            [4.2, 4.2, 5.0],   # C1
            [5.7, 4.2, 5.0],   # C2
            [7.2, 4.2, 5.0],   # C3
            [8.7, 4.2, 5.0],   # C4
            [10.2, 4.2, 5.0],  # C5
            [11.7, 4.2, 5.0],  # C6
            [13.2, 4.2, 5.0],  # C7
            [14.7, 4.2, 5.0],  # C8
            [16.2, 4.2, 5.0],  # C9 (酮基碳)
            # 酮基氧
            [17.7, 4.2, 5.0],  # O
            # H原子
            [4.2, 5.7, 5.0],   # H on C1
            [5.7, 5.7, 5.0],   # H on C2
            [7.2, 5.7, 5.0],   # H on C3
            [8.7, 5.7, 5.0],   # H on C4
            [10.2, 5.7, 5.0],  # H on C5
            [11.7, 5.7, 5.0],  # H on C6
            [13.2, 5.7, 5.0],  # H on C7
            [14.7, 5.7, 5.0],  # H on C8
            [16.2, 5.7, 5.0],  # H on C9
        ]
        
        molecule_symbols = ['C'] * 9 + ['O'] + ['H'] * 9
        
        # 合并
        all_positions = surface_positions + molecule_positions
        all_symbols = surface_symbols + molecule_symbols
        
        atoms = Atoms(symbols=all_symbols, positions=all_positions)
        
        # 设置晶胞
        atoms.set_cell([8.4, 8.4, 20.0])
        atoms.set_pbc([True, True, False])
        
        # 保存
        filename = "sample_surface_with_2_nonanone.cif"
        write(filename, atoms)
        
        print(f"   ✅ Sample structure saved as: {filename}")
        print(f"   📊 Total atoms: {len(atoms)}")
        print(f"   🏗️  Surface atoms: {len(surface_positions)}")
        print(f"   🧪 Molecule atoms: {len(molecule_positions)}")
        
        return filename
        
    except ImportError as e:
        print(f"❌ ASE not available: {e}")
        return None


def main():
    """主函数"""
    print("🚀 Surface + Molecule Structure Analysis")
    print("=" * 60)
    
    # 检查是否有结构文件
    structure_files = [
        "surface_with_molecule.cif",
        "sample_surface_with_2_nonanone.cif",
        "your_surface_structure.cif"
    ]
    
    structure_file = None
    for file in structure_files:
        if os.path.exists(file):
            structure_file = file
            break
    
    if structure_file is None:
        print("⚠️  No structure file found, creating sample...")
        structure_file = create_sample_surface_molecule()
        if not structure_file:
            print("❌ Failed to create sample structure")
            return
    
    print(f"📖 Using structure file: {structure_file}")
    
    # 分析结构
    analysis_result = analyze_surface_molecule_structure(structure_file)
    
    if analysis_result['success']:
        print(f"\n✅ Structure analysis successful!")
        print(f"   🏗️  Surface atoms: {analysis_result['surface_atoms']}")
        print(f"   🧪 Molecule atoms: {analysis_result['molecule_atoms_count']}")
        print(f"   📊 Total atoms: {analysis_result['total_atoms']}")
        print(f"   💾 Molecule file: {analysis_result['molecule_file']}")
        
        # 运行C原子优化
        try:
            response = input("\nRun carbon optimization on extracted molecule? (y/n): ").lower().strip()
            if response in ['y', 'yes', '是']:
                success = run_carbon_optimization_on_extracted_molecule(
                    analysis_result['molecule_file'],
                    temperature=400,
                    md_steps=10000,
                    opt_trials=50
                )
                
                if success:
                    print("\n🎉 Complete analysis finished!")
                    print("📁 Check 'carbon_optimization_results/' for optimization results")
                else:
                    print("\n❌ Optimization failed")
            else:
                print("\n👋 Analysis completed, optimization skipped")
                
        except KeyboardInterrupt:
            print("\n👋 User interrupted")
    else:
        print(f"\n❌ Structure analysis failed: {analysis_result['reason']}")
        print("\n💡 Tips:")
        print("  - Ensure your structure contains organic molecules (C, H, O, N, etc.)")
        print("  - Check that molecules are positioned above the surface")
        print("  - Verify the file format is supported (CIF, XYZ, etc.)")


if __name__ == "__main__":
    main()
