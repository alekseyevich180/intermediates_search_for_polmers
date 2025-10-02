#!/usr/bin/env python3
"""
修复的C原子优化分析示例脚本

这个脚本修复了导入问题，提供了多种调用方式
"""

import os
import sys
from pathlib import Path

# 添加多个可能的路径
current_dir = Path(__file__).parent
possible_paths = [
    current_dir,
    current_dir.parent / "optuna+MD",
    current_dir / ".." / "optuna+MD"
]

for path in possible_paths:
    if path.exists():
        sys.path.insert(0, str(path))

def method1_direct_command():
    """方法1: 直接使用命令行调用"""
    print("🔧 Method 1: Direct command line usage")
    print("=" * 50)
    
    # 查找search.py文件
    search_py_path = None
    for path in possible_paths:
        potential_path = path / "search.py"
        if potential_path.exists():
            search_py_path = str(potential_path.absolute())
            break
    
    if not search_py_path:
        print("❌ Cannot find search.py file!")
        return False
    
    print(f"✅ Found search.py at: {search_py_path}")
    print("\n📝 Command examples:")
    print(f"  python {search_py_path} --carbon-optimization your_structure.cif")
    print(f"  python {search_py_path} --carbon-optimization 2-nonanone.cif 400 100000 100")
    
    return True


def method2_import_with_try():
    """方法2: 尝试导入，如果失败则使用命令行"""
    print("\n🔧 Method 2: Import with fallback")
    print("=" * 50)
    
    try:
        # 尝试导入
        from search import run_carbon_optimization_analysis
        print("✅ Successfully imported from search.py")
        
        # 如果导入成功，直接使用
        print("📖 You can now use:")
        print("  result = run_carbon_optimization_analysis(")
        print("      structure_file='your_structure.cif',")
        print("      temperature=400,")
        print("      md_steps=100000,")
        print("      opt_trials=100")
        print("  )")
        
        return True
        
    except ImportError as e:
        print(f"⚠️  Import failed: {e}")
        print("🔄 Falling back to command line method...")
        return method1_direct_command()


def method3_subprocess_call():
    """方法3: 使用subprocess调用"""
    print("\n🔧 Method 3: Subprocess call")
    print("=" * 50)
    
    import subprocess
    
    # 查找search.py
    search_py_path = None
    for path in possible_paths:
        potential_path = path / "search.py"
        if potential_path.exists():
            search_py_path = str(potential_path.absolute())
            break
    
    if not search_py_path:
        print("❌ Cannot find search.py file!")
        return False
    
    # 示例调用
    structure_file = "sample_structure.cif"  # 替换为您的结构文件
    
    cmd = [
        sys.executable,
        search_py_path,
        "--carbon-optimization",
        structure_file,
        "400",    # temperature
        "10000",  # md_steps
        "50"      # opt_trials
    ]
    
    print(f"📝 Example subprocess call:")
    print(f"  cmd = {cmd}")
    print(f"  result = subprocess.run(cmd)")
    
    return True


def method4_copy_functions():
    """方法4: 复制必要的函数到当前文件"""
    print("\n🔧 Method 4: Copy functions locally")
    print("=" * 50)
    
    print("📝 You can copy the needed functions from search.py:")
    print("  1. FunctionalGroupAnalyzer")
    print("  2. MDStabilitySearcher") 
    print("  3. ProgressiveCarbonOptimizer")
    print("  4. run_carbon_optimization_analysis")
    print("\nThen use them directly in your script.")
    
    return True


def create_sample_structure():
    """创建示例结构文件"""
    print("\n🔧 Creating sample structure...")
    
    try:
        from ase import Atoms
        from ase.io import write
        
        # 简单的酮分子结构
        positions = [
            [0.0, 0.0, 0.0],    # C1
            [1.5, 0.0, 0.0],    # C2
            [3.0, 0.0, 0.0],    # C3
            [1.5, 1.5, 0.0],    # O (酮基氧)
        ]
        symbols = ['C', 'C', 'C', 'O']
        atoms = Atoms(symbols=symbols, positions=positions)
        
        # 保存结构
        write("sample_structure.cif", atoms)
        print("✅ Sample structure saved as sample_structure.cif")
        
        return True
        
    except ImportError:
        print("❌ ASE not available, cannot create sample structure")
        return False


def main():
    """主函数"""
    print("🚀 Fixed Carbon Optimization Analysis")
    print("=" * 60)
    
    print(f"📂 Current directory: {os.getcwd()}")
    print(f"📂 Script location: {Path(__file__).parent}")
    
    # 显示所有可用方法
    methods = [
        method1_direct_command,
        method2_import_with_try,
        method3_subprocess_call,
        method4_copy_functions
    ]
    
    for i, method in enumerate(methods, 1):
        try:
            method()
        except Exception as e:
            print(f"❌ Method {i} failed: {e}")
    
    # 创建示例结构
    create_sample_structure()
    
    print("\n" + "=" * 60)
    print("🎉 All methods demonstrated!")
    print("\n💡 Recommended approach:")
    print("  1. Use direct command line: python search.py --carbon-optimization file.cif")
    print("  2. Or use the standalone script: python standalone_carbon_optimization.py")


if __name__ == "__main__":
    main()
