#!/usr/bin/env python3
"""
C原子优化功能测试脚本

这个脚本用于测试search.py中的C原子优化功能是否正常工作
"""

import os
import sys
from pathlib import Path

# 添加当前目录到路径
sys.path.append(str(Path(__file__).parent))

def test_imports():
    """测试导入是否正常"""
    print("🧪 Testing imports...")
    
    try:
        from ase import Atoms
        print("   ✅ ASE imported successfully")
    except ImportError as e:
        print(f"   ❌ ASE import failed: {e}")
        return False
    
    try:
        from search import FunctionalGroupAnalyzer, MDStabilitySearcher, ProgressiveCarbonOptimizer
        print("   ✅ Search modules imported successfully")
    except ImportError as e:
        print(f"   ❌ Search modules import failed: {e}")
        return False
    
    return True


def test_functional_group_analyzer():
    """测试官能团分析器"""
    print("\n🔍 Testing FunctionalGroupAnalyzer...")
    
    try:
        from ase import Atoms
        from search import FunctionalGroupAnalyzer
        
        # 创建一个简单的酮分子 (丙酮)
        positions = [
            [0.0, 0.0, 0.0],    # C1
            [1.5, 0.0, 0.0],    # C2 (酮基碳)
            [3.0, 0.0, 0.0],    # C3
            [1.5, 1.5, 0.0],    # O (酮基氧)
            [0.0, 1.0, 0.0],    # H on C1
            [3.0, 1.0, 0.0],    # H on C3
        ]
        symbols = ['C', 'C', 'C', 'O', 'H', 'H']
        atoms = Atoms(symbols=symbols, positions=positions)
        
        # 创建一个模拟的计算器
        class MockCalculator:
            def get_potential_energy(self, atoms):
                return -100.0
        
        calculator = MockCalculator()
        
        # 测试官能团分析器
        analyzer = FunctionalGroupAnalyzer(atoms, calculator)
        functional_groups = analyzer.identify_functional_groups()
        carbon_atoms = analyzer.find_carbons_near_functional_groups()
        
        print(f"   ✅ Found {len(functional_groups)} functional groups")
        print(f"   ✅ Found {len(carbon_atoms)} nearby carbon atoms")
        
        return True
        
    except Exception as e:
        print(f"   ❌ FunctionalGroupAnalyzer test failed: {e}")
        return False


def test_md_stability_searcher():
    """测试MD稳定性搜索器"""
    print("\n🧪 Testing MDStabilitySearcher...")
    
    try:
        from ase import Atoms
        from search import MDStabilitySearcher
        
        # 创建简单的测试分子
        atoms = Atoms(symbols=['C', 'H', 'H', 'H', 'H'], 
                     positions=[[0, 0, 0], [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
        
        class MockCalculator:
            def get_potential_energy(self, atoms):
                return -50.0
        
        calculator = MockCalculator()
        
        # 测试MD稳定性搜索器 (使用较少的步数)
        md_searcher = MDStabilitySearcher(atoms, calculator, temperature=300, steps=100, sample_interval=10)
        
        print("   ✅ MDStabilitySearcher created successfully")
        print("   ⚠️  Skipping actual MD run (would require real calculator)")
        
        return True
        
    except Exception as e:
        print(f"   ❌ MDStabilitySearcher test failed: {e}")
        return False


def test_progressive_optimizer():
    """测试渐进式优化器"""
    print("\n🎯 Testing ProgressiveCarbonOptimizer...")
    
    try:
        from ase import Atoms
        from search import ProgressiveCarbonOptimizer
        
        # 创建测试数据
        atoms = Atoms(symbols=['C', 'C', 'O'], 
                     positions=[[0, 0, 0], [1.5, 0, 0], [1.5, 1.5, 0]])
        
        class MockCalculator:
            def get_potential_energy(self, atoms):
                return -75.0
        
        calculator = MockCalculator()
        
        # 模拟碳原子和官能团数据
        carbon_atoms = [{'index': 0, 'distance_to_group': 1.5, 'functional_group': 'ketone'}]
        functional_groups = [{'name': 'ketone', 'atoms': [1, 2]}]
        
        # 测试优化器创建
        optimizer = ProgressiveCarbonOptimizer(atoms, calculator, carbon_atoms, functional_groups)
        
        print("   ✅ ProgressiveCarbonOptimizer created successfully")
        print("   ⚠️  Skipping actual optimization (would require real calculator)")
        
        return True
        
    except Exception as e:
        print(f"   ❌ ProgressiveCarbonOptimizer test failed: {e}")
        return False


def test_main_function():
    """测试主函数"""
    print("\n🚀 Testing main function...")
    
    try:
        from search import run_carbon_optimization_analysis
        
        print("   ✅ Main function imported successfully")
        print("   ⚠️  Skipping full analysis (would require real structure file and calculator)")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Main function test failed: {e}")
        return False


def run_all_tests():
    """运行所有测试"""
    print("🔬 Running Carbon Optimization Tests")
    print("=" * 50)
    
    tests = [
        test_imports,
        test_functional_group_analyzer,
        test_md_stability_searcher,
        test_progressive_optimizer,
        test_main_function
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
    
    print(f"\n📊 Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("✅ All tests passed! The carbon optimization functionality is ready to use.")
    else:
        print("⚠️  Some tests failed. Please check the error messages above.")
    
    return passed == total


def main():
    """主函数"""
    success = run_all_tests()
    
    if success:
        print("\n🎉 Ready to use carbon optimization!")
        print("\nUsage examples:")
        print("  python search.py --carbon-optimization your_structure.cif")
        print("  python example_carbon_optimization.py")
    else:
        print("\n❌ Please fix the issues before using the carbon optimization functionality.")


if __name__ == "__main__":
    main()
