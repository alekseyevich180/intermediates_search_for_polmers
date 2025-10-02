#!/usr/bin/env python3
"""
氧化反应中间状态分析示例脚本

这个脚本演示如何使用改进后的md.py程序来寻找大分子在表面吸附后的氧化反应中间状态。

使用方法:
1. 准备结构文件: surface_with_molecule.cif
2. 运行: python example_oxidation.py
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
    from md import (
        analyze_oxidation_reaction, 
        find_transition_states, 
        visualize_reaction_network,
        ReactionDetector,
        IntermediateTracker
    )
    
except ImportError as e:
    print(f"❌ 导入错误: {e}")
    print("请确保已安装所需的依赖包:")
    print("  - ase")
    print("  - pfp_api_client")
    print("  - matplotlib")
    print("  - plotly")
    sys.exit(1)


def create_sample_structure():
    """创建示例结构文件（如果不存在）"""
    structure_file = "surface_with_molecule.cif"
    
    if not os.path.exists(structure_file):
        print(f"⚠️  结构文件 {structure_file} 不存在")
        print("请准备包含大分子和表面的结构文件，或修改脚本中的文件名")
        return False
    
    return True


def run_basic_oxidation_analysis():
    """运行基本的氧化反应分析"""
    print("🔥 开始氧化反应分析...")
    
    # 检查结构文件
    if not create_sample_structure():
        return False
    
    try:
        # 读取结构
        print("   📖 读取结构文件...")
        atoms = read("surface_with_molecule.cif")
        print(f"   ✅ 成功读取结构: {len(atoms)} 个原子")
        
        # 设置计算器
        print("   ⚙️  设置计算器...")
        estimator = Estimator(model_version="v8.0.0", calc_mode=EstimatorCalcMode.CRYSTAL_U0_PLUS_D3)
        calculator = ASECalculator(estimator)
        
        # 运行氧化反应分析
        print("   🧪 开始MD模拟和反应监控...")
        detector, tracker = analyze_oxidation_reaction(
            atoms, 
            calculator, 
            temperature=400,  # 较高温度促进反应
            timestep=0.5, 
            steps=5000  # 减少步数用于演示
        )
        
        # 分析结果
        print("\n📊 分析结果:")
        if detector.intermediates:
            print(f"   🎯 发现 {len(detector.intermediates)} 个反应中间体")
            
            # 显示中间体信息
            for i, intermediate in enumerate(detector.intermediates):
                print(f"     中间体 {i}: 步骤 {intermediate['step']}, 能量 {intermediate['energy']:.4f} eV")
                if intermediate['bond_changes']:
                    print(f"       键变化: {len(intermediate['bond_changes'])} 个")
            
            # 寻找过渡态
            print("\n🎯 搜索过渡态...")
            os.makedirs("transition_states", exist_ok=True)
            transition_states = find_transition_states(detector.intermediates, calculator)
            
            if transition_states:
                print(f"   ✅ 发现 {len(transition_states)} 个过渡态")
                for i, ts in enumerate(transition_states):
                    print(f"     过渡态 {i}: 能量 {ts['energy']:.4f} eV")
                
                # 创建可视化
                print("\n📈 创建反应网络可视化...")
                try:
                    fig = visualize_reaction_network(detector.intermediates, transition_states)
                    print("   ✅ 反应网络图已保存为 reaction_network.html")
                except ImportError:
                    print("   ⚠️  plotly未安装，跳过可视化")
                
            else:
                print("   ⚠️  未发现过渡态")
                
            print("\n✅ 分析完成!")
            print("   📁 结果保存在以下目录:")
            print("     - oxidation_analysis/ (中间体结构)")
            print("     - transition_states/ (过渡态结构)")
            print("     - reaction_network.html (可视化)")
            
        else:
            print("   ⚠️  未检测到反应中间体")
            print("   建议:")
            print("     - 提高模拟温度")
            print("     - 增加模拟步数")
            print("     - 检查分子是否包含可反应的基团")
        
        return True
        
    except Exception as e:
        print(f"❌ 分析过程中出现错误: {e}")
        return False


def run_detailed_analysis():
    """运行详细分析（需要更多计算资源）"""
    print("\n🔬 运行详细分析...")
    
    if not create_sample_structure():
        return False
    
    try:
        # 读取结构
        atoms = read("surface_with_molecule.cif")
        
        # 设置计算器
        estimator = Estimator(model_version="v8.0.0", calc_mode=EstimatorCalcMode.CRYSTAL_U0_PLUS_D3)
        calculator = ASECalculator(estimator)
        
        # 创建自定义监控器
        detector = ReactionDetector(atoms, calculator, bond_threshold=1.8)  # 更严格的键长阈值
        tracker = IntermediateTracker("detailed_analysis")
        
        print("   ⚙️  使用更严格的参数进行分析...")
        print("   📊 键长阈值: 1.8 Å")
        print("   🔍 将进行更详细的键变化分析")
        
        # 这里可以添加更详细的分析逻辑
        # 例如：分析特定的反应类型、键角变化等
        
        return True
        
    except Exception as e:
        print(f"❌ 详细分析失败: {e}")
        return False


def main():
    """主函数"""
    print("🚀 氧化反应中间状态分析示例")
    print("=" * 50)
    
    # 检查工作目录
    print(f"📂 工作目录: {os.getcwd()}")
    
    # 运行基本分析
    success = run_basic_oxidation_analysis()
    
    if success:
        print("\n" + "=" * 50)
        print("🎉 基本分析完成!")
        
        # 询问是否运行详细分析
        try:
            response = input("\n是否运行详细分析? (y/n): ").lower().strip()
            if response in ['y', 'yes', '是']:
                run_detailed_analysis()
        except KeyboardInterrupt:
            print("\n👋 用户取消操作")
    
    else:
        print("\n❌ 分析失败")
        print("\n💡 故障排除建议:")
        print("1. 检查结构文件是否存在且格式正确")
        print("2. 确认计算器配置正确")
        print("3. 检查依赖包是否完整安装")
        print("4. 查看错误信息进行调试")


if __name__ == "__main__":
    main()
