#!/usr/bin/env python3
"""
C原子优化分析示例脚本

这个脚本演示如何使用search.py程序来：
1. 识别官能团周围的C原子
2. 通过MD模拟找到最稳定结构
3. 使用Optuna从近到远优化C原子位置

使用方法:
1. 准备反应中间体结构文件 (例如: 2-nonanone.cif)
2. 运行: python example_carbon_optimization.py
"""

import os
import sys
from pathlib import Path

# 添加当前目录到路径
sys.path.append(str(Path(__file__).parent))

try:
    from ase.io import read, write
    from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
    from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
    
    # 导入我们的分析函数
    from search import (
        FunctionalGroupAnalyzer,
        MDStabilitySearcher,
        ProgressiveCarbonOptimizer,
        run_carbon_optimization_analysis
    )
    
except ImportError as e:
    print(f"❌ 导入错误: {e}")
    print("请确保已安装所需的依赖包:")
    print("  - ase")
    print("  - pfp_api_client")
    print("  - optuna")
    print("  - scipy")
    print("  - sklearn")
    sys.exit(1)


def create_sample_2_nonanone():
    """创建一个示例的2-nonanone分子结构"""
    print("🔧 Creating sample 2-nonanone structure...")
    
    # 简化的2-nonanone结构 (C9H18O)
    # 这是一个示例，实际使用时请提供真实的结构文件
    
    from ase import Atoms
    
    # 2-nonanone的基本结构
    positions = [
        # 主链碳原子 (简化的线性结构)
        [0.0, 0.0, 0.0],      # C1
        [1.5, 0.0, 0.0],      # C2
        [3.0, 0.0, 0.0],      # C3
        [4.5, 0.0, 0.0],      # C4
        [6.0, 0.0, 0.0],      # C5
        [7.5, 0.0, 0.0],      # C6
        [9.0, 0.0, 0.0],      # C7
        [10.5, 0.0, 0.0],     # C8
        [12.0, 0.0, 0.0],     # C9 (酮基碳)
        # 酮基氧原子
        [13.5, 0.0, 0.0],     # O (酮基氧)
        # 氢原子 (简化的位置)
        [0.0, 1.5, 0.0],      # H on C1
        [1.5, 1.5, 0.0],      # H on C2
        [3.0, 1.5, 0.0],      # H on C3
        [4.5, 1.5, 0.0],      # H on C4
        [6.0, 1.5, 0.0],      # H on C5
        [7.5, 1.5, 0.0],      # H on C6
        [9.0, 1.5, 0.0],      # H on C7
        [10.5, 1.5, 0.0],     # H on C8
        [12.0, 1.5, 0.0],     # H on C9
    ]
    
    symbols = ['C'] * 9 + ['O'] + ['H'] * 9
    
    atoms = Atoms(symbols=symbols, positions=positions)
    
    # 保存示例结构
    write("sample_2_nonanone.cif", atoms)
    print("   Sample structure saved as sample_2_nonanone.cif")
    
    return atoms


def run_basic_carbon_analysis():
    """运行基本的C原子分析"""
    print("🔬 Starting basic carbon analysis...")
    
    # 检查是否有结构文件，如果没有则创建示例
    structure_files = ["2-nonanone.cif", "sample_2_nonanone.cif"]
    structure_file = None
    
    for file in structure_files:
        if os.path.exists(file):
            structure_file = file
            break
    
    if structure_file is None:
        print("⚠️  No structure file found, creating sample structure...")
        atoms = create_sample_2_nonanone()
        structure_file = "sample_2_nonanone.cif"
    
    print(f"📖 Using structure file: {structure_file}")
    
    try:
        # 运行C原子优化分析
        result = run_carbon_optimization_analysis(
            structure_file=structure_file,
            temperature=400,  # 较高温度促进结构探索
            md_steps=10000,   # 减少步数用于演示
            opt_trials=50     # 减少试验次数用于演示
        )
        
        if result:
            print("\n✅ Analysis completed successfully!")
            
            # 显示详细结果
            summary = result['summary']
            print(f"\n📊 Detailed Results:")
            print(f"   Energy improvement: {summary['energy_improvement']:.4f} eV")
            print(f"   Functional groups found: {len(summary['functional_groups'])}")
            print(f"   Carbon atoms optimized: {len(summary['carbon_atoms'])}")
            
            # 显示官能团信息
            if summary['functional_groups']:
                print(f"\n🔍 Functional Groups:")
                for group in summary['functional_groups']:
                    print(f"   - {group['name']}: atoms {group['atoms']}")
            
            # 显示C原子信息
            if summary['carbon_atoms']:
                print(f"\n⚛️  Carbon Atoms (first 5):")
                for i, carbon in enumerate(summary['carbon_atoms'][:5]):
                    print(f"   - Carbon {carbon['index']}: distance {carbon['distance_to_group']:.2f} Å to {carbon['functional_group']}")
            
            print(f"\n📁 Results saved in: carbon_optimization_results/")
            
        else:
            print("❌ Analysis failed")
            
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()


def run_detailed_analysis():
    """运行详细分析（需要更多计算资源）"""
    print("\n🔬 Running detailed analysis...")
    
    structure_file = input("Enter structure file path (or press Enter for sample): ").strip()
    if not structure_file:
        structure_file = "sample_2_nonanone.cif"
    
    if not os.path.exists(structure_file):
        print(f"❌ File not found: {structure_file}")
        return
    
    try:
        # 获取用户参数
        print("\n⚙️  Configuration:")
        temperature = input("MD temperature (K) [400]: ").strip()
        temperature = int(temperature) if temperature else 400
        
        md_steps = input("MD steps [100000]: ").strip()
        md_steps = int(md_steps) if md_steps else 100000
        
        opt_trials = input("Optimization trials per carbon [100]: ").strip()
        opt_trials = int(opt_trials) if opt_trials else 100
        
        print(f"\n🚀 Starting analysis with:")
        print(f"   Temperature: {temperature} K")
        print(f"   MD steps: {md_steps}")
        print(f"   Optimization trials: {opt_trials}")
        
        # 运行分析
        result = run_carbon_optimization_analysis(
            structure_file=structure_file,
            temperature=temperature,
            md_steps=md_steps,
            opt_trials=opt_trials
        )
        
        if result:
            print("\n✅ Detailed analysis completed!")
            print("   Check carbon_optimization_results/ for detailed results")
        else:
            print("❌ Detailed analysis failed")
            
    except KeyboardInterrupt:
        print("\n👋 Analysis interrupted by user")
    except Exception as e:
        print(f"❌ Error during detailed analysis: {e}")


def demonstrate_individual_components():
    """演示各个组件的独立使用"""
    print("\n🧪 Demonstrating individual components...")
    
    try:
        # 读取结构
        structure_file = "sample_2_nonanone.cif"
        if not os.path.exists(structure_file):
            create_sample_2_nonanone()
        
        atoms = read(structure_file)
        print(f"📖 Loaded structure: {len(atoms)} atoms")
        
        # 设置计算器
        estimator = Estimator(model_version="v8.0.0", calc_mode=EstimatorCalcMode.CRYSTAL_U0_PLUS_D3)
        calculator = ASECalculator(estimator)
        
        # 1. 官能团分析
        print("\n🔍 Functional Group Analysis:")
        analyzer = FunctionalGroupAnalyzer(atoms, calculator)
        functional_groups = analyzer.identify_functional_groups()
        print(f"   Found {len(functional_groups)} functional groups")
        
        carbon_atoms = analyzer.find_carbons_near_functional_groups()
        print(f"   Found {len(carbon_atoms)} nearby carbon atoms")
        
        # 2. MD稳定性搜索
        print("\n🧪 MD Stability Search:")
        md_searcher = MDStabilitySearcher(atoms, calculator, temperature=300, steps=1000, sample_interval=100)
        most_stable = md_searcher.find_most_stable_structure()
        print(f"   Most stable energy: {most_stable['energy']:.4f} eV")
        
        # 3. 渐进式优化
        if carbon_atoms:
            print("\n🎯 Progressive Optimization:")
            optimizer = ProgressiveCarbonOptimizer(atoms, calculator, carbon_atoms[:3], functional_groups)  # 只优化前3个C原子
            history = optimizer.optimize_carbons_progressively(n_trials=20)  # 减少试验次数
            print(f"   Optimized {len(history)} carbon atoms")
        
        print("\n✅ Component demonstration completed!")
        
    except Exception as e:
        print(f"❌ Error during component demonstration: {e}")


def main():
    """主函数"""
    print("🚀 Carbon Optimization Analysis Example")
    print("=" * 50)
    
    # 检查工作目录
    print(f"📂 Working directory: {os.getcwd()}")
    
    # 运行基本分析
    run_basic_carbon_analysis()
    
    print("\n" + "=" * 50)
    print("🎉 Basic analysis completed!")
    
    # 询问是否运行详细分析
    try:
        response = input("\nRun detailed analysis? (y/n): ").lower().strip()
        if response in ['y', 'yes', '是']:
            run_detailed_analysis()
        
        # 询问是否演示各个组件
        response = input("\nDemonstrate individual components? (y/n): ").lower().strip()
        if response in ['y', 'yes', '是']:
            demonstrate_individual_components()
            
    except KeyboardInterrupt:
        print("\n👋 User interrupted")


if __name__ == "__main__":
    main()
