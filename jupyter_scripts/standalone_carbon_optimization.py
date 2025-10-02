#!/usr/bin/env python3
"""
独立的C原子优化分析脚本

这个脚本不依赖导入search.py，直接调用命令行接口
使用方法:
1. 准备反应中间体结构文件 (例如: 2-nonanone.cif)
2. 运行: python standalone_carbon_optimization.py
"""

import os
import sys
import subprocess
from pathlib import Path

def create_sample_surface_with_molecule():
    """创建一个包含表面和分子的示例结构"""
    print("🔧 Creating sample surface with molecule structure...")
    
    try:
        from ase import Atoms
        from ase.io import write
        import numpy as np
        
        # 创建表面 (简化的金属表面)
        surface_positions = []
        surface_symbols = []
        
        # 创建一个3x3的金属表面
        for i in range(3):
            for j in range(3):
                surface_positions.append([i*2.5, j*2.5, 0.0])  # 底层
                surface_positions.append([i*2.5, j*2.5, 2.0])  # 顶层
                surface_symbols.extend(['Pt', 'Pt'])  # 假设是铂表面
        
        # 创建分子 (2-nonanone)
        molecule_positions = [
            # 主链碳原子 (在表面上方的z位置)
            [3.0, 3.0, 4.0],      # C1
            [4.5, 3.0, 4.0],      # C2
            [6.0, 3.0, 4.0],      # C3
            [7.5, 3.0, 4.0],      # C4
            [9.0, 3.0, 4.0],      # C5
            [10.5, 3.0, 4.0],     # C6
            [12.0, 3.0, 4.0],     # C7
            [13.5, 3.0, 4.0],     # C8
            [15.0, 3.0, 4.0],     # C9 (酮基碳)
            # 酮基氧原子
            [16.5, 3.0, 4.0],     # O (酮基氧)
            # 氢原子
            [3.0, 4.5, 4.0],      # H on C1
            [4.5, 4.5, 4.0],      # H on C2
            [6.0, 4.5, 4.0],      # H on C3
            [7.5, 4.5, 4.0],      # H on C4
            [9.0, 4.5, 4.0],      # H on C5
            [10.5, 4.5, 4.0],     # H on C6
            [12.0, 4.5, 4.0],     # H on C7
            [13.5, 4.5, 4.0],     # H on C8
            [15.0, 4.5, 4.0],     # H on C9
        ]
        
        molecule_symbols = ['C'] * 9 + ['O'] + ['H'] * 9
        
        # 合并表面和分子
        all_positions = surface_positions + molecule_positions
        all_symbols = surface_symbols + molecule_symbols
        
        atoms = Atoms(symbols=all_symbols, positions=all_positions)
        
        # 设置周期性边界条件 (表面)
        atoms.set_cell([7.5, 7.5, 20.0])  # 足够大的z方向容纳分子
        atoms.set_pbc([True, True, False])  # x,y方向周期性，z方向非周期性
        
        # 保存示例结构
        write("sample_surface_with_molecule.cif", atoms)
        print("   Sample surface+molecule structure saved as sample_surface_with_molecule.cif")
        
        return "sample_surface_with_molecule.cif"
        
    except ImportError as e:
        print(f"❌ ASE not available: {e}")
        return None


def analyze_structure_with_surface(structure_file):
    """分析包含表面的结构文件，识别分子部分"""
    print(f"🔍 Analyzing structure with surface: {structure_file}")
    
    try:
        from ase.io import read
        from ase import Atoms
        import numpy as np
        
        # 读取完整结构
        full_atoms = read(structure_file)
        print(f"   Total atoms in structure: {len(full_atoms)}")
        
        # 分析原子类型
        symbols = full_atoms.get_chemical_symbols()
        unique_symbols = list(set(symbols))
        print(f"   Element types: {unique_symbols}")
        
        # 统计各元素数量
        element_counts = {}
        for symbol in unique_symbols:
            element_counts[symbol] = symbols.count(symbol)
        print(f"   Element counts: {element_counts}")
        
        # 识别可能的表面原子 (通常密度较高，位置较低)
        positions = full_atoms.get_positions()
        
        # 按z坐标分层分析
        z_coords = positions[:, 2]
        z_min, z_max = z_coords.min(), z_coords.max()
        z_range = z_max - z_min
        
        print(f"   Z coordinate range: {z_min:.2f} to {z_max:.2f} Å")
        
        # 识别表面层 (通常在最底部)
        surface_threshold = z_min + z_range * 0.3  # 底部30%认为是表面
        surface_indices = np.where(z_coords < surface_threshold)[0]
        
        # 识别分子层 (在表面上方的部分)
        molecule_threshold = z_min + z_range * 0.4  # 40%以上认为是分子
        molecule_indices = np.where(z_coords > molecule_threshold)[0]
        
        print(f"   Surface atoms (z < {surface_threshold:.2f}): {len(surface_indices)}")
        print(f"   Molecule atoms (z > {molecule_threshold:.2f}): {len(molecule_indices)}")
        
        # 分析分子部分的元素组成
        if len(molecule_indices) > 0:
            molecule_symbols = [symbols[i] for i in molecule_indices]
            molecule_elements = list(set(molecule_symbols))
            molecule_counts = {}
            for element in molecule_elements:
                molecule_counts[element] = molecule_symbols.count(element)
            
            print(f"   Molecule elements: {molecule_elements}")
            print(f"   Molecule element counts: {molecule_counts}")
            
            # 检查是否包含有机分子特征 (C, H, O, N等)
            organic_elements = ['C', 'H', 'O', 'N', 'S', 'P']
            has_organic = any(elem in molecule_elements for elem in organic_elements)
            
            if has_organic:
                print("   ✅ Organic molecule detected in structure")
                
                # 提取分子部分
                molecule_positions = positions[molecule_indices]
                molecule_atoms = Atoms(symbols=molecule_symbols, positions=molecule_positions)
                
                # 保存分子部分
                molecule_file = structure_file.replace('.cif', '_molecule_only.cif')
                write(molecule_file, molecule_atoms)
                print(f"   💾 Molecule-only structure saved as: {molecule_file}")
                
                return molecule_file, molecule_atoms, molecule_indices
            else:
                print("   ⚠️  No organic molecule detected")
                return None, None, None
        else:
            print("   ⚠️  No molecule layer detected")
            return None, None, None
            
    except Exception as e:
        print(f"❌ Error analyzing structure: {e}")
        return None, None, None


def find_search_py():
    """查找search.py文件的位置"""
    current_dir = Path(__file__).parent
    possible_paths = [
        current_dir / "search.py",
        current_dir.parent / "optuna+MD" / "search.py",
        current_dir / ".." / "optuna+MD" / "search.py"
    ]
    
    for path in possible_paths:
        if path.exists():
            return str(path.absolute())
    
    return None


def run_carbon_optimization_command(structure_file, temperature=400, md_steps=10000, opt_trials=50):
    """运行C原子优化命令"""
    print(f"🚀 Running carbon optimization for: {structure_file}")
    
    # 查找search.py文件
    search_py_path = find_search_py()
    if not search_py_path:
        print("❌ Cannot find search.py file!")
        print("Please ensure search.py is in one of these locations:")
        print("  - Current directory")
        print("  - ../optuna+MD/")
        return False
    
    print(f"📂 Found search.py at: {search_py_path}")
    
    # 构建命令
    cmd = [
        sys.executable,  # Python解释器
        search_py_path,
        "--carbon-optimization",
        structure_file,
        str(temperature),
        str(md_steps),
        str(opt_trials)
    ]
    
    print(f"🔧 Command: {' '.join(cmd)}")
    
    try:
        # 运行命令
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)  # 1小时超时
        
        # 显示输出
        if result.stdout:
            print("📊 Output:")
            print(result.stdout)
        
        if result.stderr:
            print("⚠️  Errors/Warnings:")
            print(result.stderr)
        
        if result.returncode == 0:
            print("✅ Carbon optimization completed successfully!")
            return True
        else:
            print(f"❌ Command failed with return code: {result.returncode}")
            return False
            
    except subprocess.TimeoutExpired:
        print("⏰ Command timed out after 1 hour")
        return False
    except Exception as e:
        print(f"❌ Error running command: {e}")
        return False


def run_basic_analysis():
    """运行基本分析"""
    print("🔬 Starting basic carbon analysis...")
    
    # 检查是否有结构文件（包括表面+分子结构）
    structure_files = [
        "surface_with_molecule.cif",
        "sample_surface_with_molecule.cif", 
        "2-nonanone.cif", 
        "sample_2_nonanone.cif"
    ]
    structure_file = None
    
    for file in structure_files:
        if os.path.exists(file):
            structure_file = file
            break
    
    if structure_file is None:
        print("⚠️  No structure file found, creating sample surface+molecule structure...")
        structure_file = create_sample_surface_with_molecule()
        if not structure_file:
            return False
    
    print(f"📖 Using structure file: {structure_file}")
    
    # 分析结构，识别分子部分
    molecule_file, molecule_atoms, molecule_indices = analyze_structure_with_surface(structure_file)
    
    if molecule_file:
        print(f"✅ Molecule extracted successfully!")
        print(f"   Molecule file: {molecule_file}")
        print(f"   Molecule atoms: {len(molecule_atoms)}")
        
        # 使用提取的分子部分进行优化
        target_file = molecule_file
    else:
        print("⚠️  Using full structure (assuming it's molecule-only)")
        target_file = structure_file
    
    # 运行优化（使用较小的参数进行演示）
    success = run_carbon_optimization_command(
        structure_file=target_file,
        temperature=400,
        md_steps=5000,  # 减少步数用于演示
        opt_trials=30   # 减少试验次数用于演示
    )
    
    if success:
        print("\n📁 Check 'carbon_optimization_results/' directory for results")
        print("   - final_optimized_structure.cif")
        print("   - optimization_summary.json")
    
    return success


def run_detailed_analysis():
    """运行详细分析"""
    print("\n🔬 Running detailed analysis...")
    
    # 获取用户输入
    structure_file = input("Enter structure file path (or press Enter for sample): ").strip()
    if not structure_file:
        structure_file = "sample_2_nonanone.cif"
        if not os.path.exists(structure_file):
            structure_file = create_sample_2_nonanone()
            if not structure_file:
                return False
    
    if not os.path.exists(structure_file):
        print(f"❌ File not found: {structure_file}")
        return False
    
    try:
        # 获取参数
        print("\n⚙️  Configuration:")
        temperature = input("MD temperature (K) [400]: ").strip()
        temperature = int(temperature) if temperature else 400
        
        md_steps = input("MD steps [100000]: ").strip()
        md_steps = int(md_steps) if md_steps else 100000
        
        opt_trials = input("Optimization trials per carbon [100]: ").strip()
        opt_trials = int(opt_trials) if opt_trials else 100
        
        print(f"\n🚀 Starting detailed analysis:")
        print(f"   Structure: {structure_file}")
        print(f"   Temperature: {temperature} K")
        print(f"   MD steps: {md_steps}")
        print(f"   Optimization trials: {opt_trials}")
        
        # 运行分析
        success = run_carbon_optimization_command(
            structure_file=structure_file,
            temperature=temperature,
            md_steps=md_steps,
            opt_trials=opt_trials
        )
        
        if success:
            print("\n✅ Detailed analysis completed!")
            print("   Check carbon_optimization_results/ for detailed results")
        
        return success
        
    except KeyboardInterrupt:
        print("\n👋 Analysis interrupted by user")
        return False
    except Exception as e:
        print(f"❌ Error during detailed analysis: {e}")
        return False


def show_usage_examples():
    """显示使用示例"""
    print("\n📚 Usage Examples:")
    print("=" * 50)
    
    search_py_path = find_search_py()
    if search_py_path:
        print(f"Python script found at: {search_py_path}")
        print("\nCommand line examples:")
        print(f"  python {search_py_path} --carbon-optimization your_structure.cif")
        print(f"  python {search_py_path} --carbon-optimization 2-nonanone.cif 400 100000 100")
        print(f"  python {search_py_path} --carbon-optimization sample_2_nonanone.cif 300 50000 75")
    else:
        print("❌ search.py not found!")
    
    print("\nParameter explanation:")
    print("  --carbon-optimization <file> [temperature] [md_steps] [opt_trials]")
    print("    file: Structure file path (required)")
    print("    temperature: MD simulation temperature in K (default: 300)")
    print("    md_steps: Number of MD steps (default: 100000)")
    print("    opt_trials: Optimization trials per carbon (default: 100)")


def main():
    """主函数"""
    print("🚀 Standalone Carbon Optimization Analysis")
    print("=" * 50)
    
    # 检查工作目录
    print(f"📂 Working directory: {os.getcwd()}")
    
    # 检查search.py是否存在
    search_py_path = find_search_py()
    if not search_py_path:
        print("❌ Cannot find search.py file!")
        print("Please ensure search.py is available.")
        show_usage_examples()
        return
    
    print(f"✅ Found search.py at: {search_py_path}")
    
    # 运行基本分析
    success = run_basic_analysis()
    
    if success:
        print("\n" + "=" * 50)
        print("🎉 Basic analysis completed!")
        
        # 询问是否运行详细分析
        try:
            response = input("\nRun detailed analysis? (y/n): ").lower().strip()
            if response in ['y', 'yes', '是']:
                run_detailed_analysis()
            
            # 显示使用示例
            response = input("\nShow usage examples? (y/n): ").lower().strip()
            if response in ['y', 'yes', '是']:
                show_usage_examples()
                
        except KeyboardInterrupt:
            print("\n👋 User interrupted")
    else:
        print("\n❌ Basic analysis failed")
        show_usage_examples()


if __name__ == "__main__":
    main()
