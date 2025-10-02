#!/usr/bin/env python3
"""
分析主控制脚本

功能: 统一控制所有分析模块，提供完整的分析流程
输入: 结构文件 (支持表面+分子或纯分子)
输出: 完整的分析报告和结果文件
"""

import argparse
import os
import sys
import subprocess
import time
from pathlib import Path

def run_analysis_module(module_name, structure_file, output_dir, **kwargs):
    """运行分析模块"""
    print(f"\n{'='*60}")
    print(f"🚀 Running {module_name}")
    print(f"{'='*60}")
    
    # 构建命令
    cmd = [sys.executable, f"{module_name}.py", structure_file]
    
    # 添加输出目录参数
    if module_name != "surface_extractor":
        cmd.extend(["--output_dir", output_dir])
    
    # 添加其他参数
    for key, value in kwargs.items():
        cmd.extend([f"--{key}", str(value)])
    
    print(f"🔧 Command: {' '.join(cmd)}")
    
    try:
        # 运行模块
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        end_time = time.time()
        
        print(f"⏱️  Execution time: {end_time - start_time:.2f} seconds")
        
        if result.stdout:
            print("📊 Output:")
            print(result.stdout)
        
        if result.stderr:
            print("⚠️  Warnings/Errors:")
            print(result.stderr)
        
        if result.returncode == 0:
            print(f"✅ {module_name} completed successfully!")
            return True
        else:
            print(f"❌ {module_name} failed with return code: {result.returncode}")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"⏰ {module_name} timed out after 1 hour")
        return False
    except Exception as e:
        print(f"❌ Error running {module_name}: {e}")
        return False


def run_complete_analysis(structure_file, analysis_modules=None, output_dir="complete_analysis_results"):
    """运行完整的分析流程"""
    print("🔬 Complete Structure Analysis Pipeline")
    print("=" * 80)
    print(f"📖 Structure file: {structure_file}")
    print(f"📁 Output directory: {output_dir}")
    
    if analysis_modules is None:
        analysis_modules = [
            "surface_extractor",
            "functional_group_analyzer", 
            "md_stability_searcher",
            "carbon_optimizer",
            "reaction_detector"
        ]
    
    print(f"🔧 Analysis modules: {', '.join(analysis_modules)}")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 存储结果
    results = {
        'structure_file': structure_file,
        'analysis_modules': analysis_modules,
        'output_dir': output_dir,
        'start_time': time.time(),
        'module_results': {},
        'successful_modules': [],
        'failed_modules': []
    }
    
    # 运行各个模块
    for module in analysis_modules:
        print(f"\n🎯 Processing module: {module}")
        
        # 为每个模块创建子目录
        module_output_dir = f"{output_dir}/{module}_results"
        
        # 设置模块特定参数
        module_kwargs = {}
        if module == "md_stability_searcher":
            module_kwargs = {
                "temperature": 300,
                "steps": 10000,
                "sample_interval": 100
            }
        elif module == "carbon_optimizer":
            module_kwargs = {
                "opt_trials": 50
            }
        elif module == "functional_group_analyzer":
            module_kwargs = {
                "max_distance": 3.0
            }
        
        # 运行模块
        success = run_analysis_module(
            module, 
            structure_file, 
            module_output_dir,
            **module_kwargs
        )
        
        # 记录结果
        results['module_results'][module] = {
            'success': success,
            'output_dir': module_output_dir,
            'parameters': module_kwargs
        }
        
        if success:
            results['successful_modules'].append(module)
        else:
            results['failed_modules'].append(module)
    
    # 完成时间
    results['end_time'] = time.time()
    results['total_time'] = results['end_time'] - results['start_time']
    
    # 生成总结报告
    generate_summary_report(results)
    
    return results


def generate_summary_report(results):
    """生成总结报告"""
    print(f"\n{'='*80}")
    print(f"📊 ANALYSIS SUMMARY REPORT")
    print(f"{'='*80}")
    
    print(f"📖 Structure file: {results['structure_file']}")
    print(f"⏱️  Total analysis time: {results['total_time']:.2f} seconds")
    print(f"📁 Output directory: {results['output_dir']}")
    
    print(f"\n✅ Successful modules ({len(results['successful_modules'])}):")
    for module in results['successful_modules']:
        print(f"   - {module}")
    
    if results['failed_modules']:
        print(f"\n❌ Failed modules ({len(results['failed_modules'])}):")
        for module in results['failed_modules']:
            print(f"   - {module}")
    
    print(f"\n📁 Module output directories:")
    for module, module_result in results['module_results'].items():
        status = "✅" if module_result['success'] else "❌"
        print(f"   {status} {module}: {module_result['output_dir']}")
    
    # 保存总结报告
    import json
    summary_file = f"{results['output_dir']}/analysis_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n💾 Summary report saved: {summary_file}")


def create_sample_structure():
    """创建示例结构文件"""
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
    parser = argparse.ArgumentParser(description="Complete Structure Analysis Pipeline")
    parser.add_argument("structure_file", nargs='?', help="Structure file path (CIF, XYZ, etc.)")
    parser.add_argument("--output_dir", type=str, default="complete_analysis_results",
                       help="Output directory name (default: complete_analysis_results)")
    parser.add_argument("--modules", nargs='+', 
                       choices=["surface_extractor", "functional_group_analyzer", 
                               "md_stability_searcher", "carbon_optimizer", "reaction_detector"],
                       help="Specific modules to run (default: all)")
    parser.add_argument("--create_sample", action="store_true",
                       help="Create a sample structure file")
    
    args = parser.parse_args()
    
    # 创建示例结构
    if args.create_sample:
        sample_file = create_sample_structure()
        if sample_file:
            print(f"\n✅ Sample structure created: {sample_file}")
            print("You can now run the analysis on this sample structure.")
        return
    
    # 检查结构文件
    if not args.structure_file:
        print("❌ Please provide a structure file or use --create_sample to create one")
        return
    
    if not os.path.exists(args.structure_file):
        print(f"❌ Structure file not found: {args.structure_file}")
        return
    
    # 检查模块文件是否存在
    available_modules = []
    if args.modules:
        modules_to_run = args.modules
    else:
        modules_to_run = ["surface_extractor", "functional_group_analyzer", 
                         "md_stability_searcher", "carbon_optimizer", "reaction_detector"]
    
    for module in modules_to_run:
        module_file = f"{module}.py"
        if os.path.exists(module_file):
            available_modules.append(module)
        else:
            print(f"⚠️  Module file not found: {module_file}")
    
    if not available_modules:
        print("❌ No analysis modules found!")
        return
    
    print(f"🔧 Available modules: {', '.join(available_modules)}")
    
    # 运行完整分析
    results = run_complete_analysis(
        args.structure_file,
        analysis_modules=available_modules,
        output_dir=args.output_dir
    )
    
    # 最终总结
    print(f"\n{'='*80}")
    if results['successful_modules']:
        print(f"🎉 Analysis completed successfully!")
        print(f"✅ {len(results['successful_modules'])}/{len(available_modules)} modules succeeded")
    else:
        print(f"❌ Analysis failed!")
        print(f"❌ {len(results['failed_modules'])}/{len(available_modules)} modules failed")
    
    print(f"📁 Results saved in: {results['output_dir']}/")
    print(f"⏱️  Total time: {results['total_time']:.2f} seconds")


if __name__ == "__main__":
    main()
