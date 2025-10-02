#!/usr/bin/env python3
"""
MD稳定性搜索独立脚本

功能: 通过MD模拟寻找最稳定的结构
输入: 结构文件 (支持表面+分子或纯分子)
输出: 最稳定结构、采样轨迹、能量分析
"""

import argparse
import os
import sys
import time
from pathlib import Path

try:
    from ase.io import read, write
    from ase import Atoms, units
    from ase.md.npt import NPT
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
    from ase.optimize import LBFGS
    from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
    from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
    import numpy as np
except ImportError as e:
    print(f"❌ 导入错误: {e}")
    print("请确保已安装所需的依赖包")
    sys.exit(1)


class MDStabilitySearcher:
    """通过MD模拟寻找最稳定的结构"""
    
    def __init__(self, atoms, calculator, temperature=300, steps=100000, sample_interval=100):
        self.atoms = atoms
        self.calculator = calculator
        self.temperature = temperature
        self.steps = steps
        self.sample_interval = sample_interval
        self.sampled_structures = []
        self.energies = []
        
    def run_md_sampling(self):
        """运行MD模拟并采样结构"""
        print(f"🧪 Running MD sampling: {self.steps} steps at {self.temperature}K")
        
        # 设置计算器
        self.atoms.calc = self.calculator
        
        # 设置初始速度
        MaxwellBoltzmannDistribution(self.atoms, temperature_K=self.temperature, force_temp=True)
        Stationary(self.atoms)
        
        # 创建MD动力学
        dyn = NPT(
            self.atoms,
            0.5 * units.fs,
            temperature_K=self.temperature,
            externalstress=0,
            ttime=100 * units.fs,
            pfactor=None,
            trajectory="md_sampling.xyz",
            loginterval=self.sample_interval,
        )
        
        # 采样回调函数
        def sample_structure():
            step = dyn.get_number_of_steps()
            if step % self.sample_interval == 0:
                energy = self.atoms.get_potential_energy()
                structure = self.atoms.copy()
                
                self.sampled_structures.append({
                    'step': step,
                    'atoms': structure,
                    'energy': energy,
                    'positions': structure.get_positions().copy()
                })
                self.energies.append(energy)
                
                if step % (self.sample_interval * 10) == 0:
                    print(f"   Sampled at step {step}, energy: {energy:.4f} eV")
        
        # 附加采样
        dyn.attach(sample_structure, interval=self.sample_interval)
        
        # 运行MD
        start_time = time.time()
        dyn.run(self.steps)
        end_time = time.time()
        
        print(f"   Sampling complete: {len(self.sampled_structures)} structures sampled")
        print(f"   MD simulation time: {end_time - start_time:.2f} seconds")
        
        return self.sampled_structures
    
    def find_most_stable_structure(self):
        """找到最稳定的结构"""
        if not self.sampled_structures:
            self.run_md_sampling()
        
        # 找到能量最低的结构
        min_energy_idx = np.argmin(self.energies)
        most_stable = self.sampled_structures[min_energy_idx]
        
        print(f"🎯 Most stable structure found at step {most_stable['step']}")
        print(f"   Energy: {most_stable['energy']:.4f} eV")
        
        return most_stable
    
    def analyze_energy_landscape(self):
        """分析能量景观"""
        if not self.energies:
            return None
        
        energies = np.array(self.energies)
        
        analysis = {
            'min_energy': float(np.min(energies)),
            'max_energy': float(np.max(energies)),
            'mean_energy': float(np.mean(energies)),
            'std_energy': float(np.std(energies)),
            'energy_range': float(np.max(energies) - np.min(energies)),
            'total_samples': len(energies)
        }
        
        return analysis


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


def run_md_stability_search(structure_file, temperature=300, steps=100000, sample_interval=100):
    """运行MD稳定性搜索"""
    print(f"🧪 Running MD stability search for: {structure_file}")
    
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
        
        # 创建MD稳定性搜索器
        searcher = MDStabilitySearcher(
            atoms, 
            calculator, 
            temperature=temperature, 
            steps=steps, 
            sample_interval=sample_interval
        )
        
        # 运行MD采样
        print(f"\n🧪 Running MD simulation...")
        sampled_structures = searcher.run_md_sampling()
        
        # 找到最稳定结构
        print(f"\n🎯 Finding most stable structure...")
        most_stable = searcher.find_most_stable_structure()
        
        # 分析能量景观
        print(f"\n📊 Analyzing energy landscape...")
        energy_analysis = searcher.analyze_energy_landscape()
        
        if energy_analysis:
            print(f"   Min energy: {energy_analysis['min_energy']:.4f} eV")
            print(f"   Max energy: {energy_analysis['max_energy']:.4f} eV")
            print(f"   Mean energy: {energy_analysis['mean_energy']:.4f} eV")
            print(f"   Energy range: {energy_analysis['energy_range']:.4f} eV")
            print(f"   Total samples: {energy_analysis['total_samples']}")
        
        # 保存结果
        output_dir = "md_stability_results"
        os.makedirs(output_dir, exist_ok=True)
        
        # 保存最稳定结构
        write(f"{output_dir}/most_stable_structure.cif", most_stable['atoms'])
        
        # 保存初始结构
        write(f"{output_dir}/initial_structure.cif", atoms)
        
        # 保存分析报告
        analysis_report = {
            'structure_file': structure_file,
            'initial_energy': initial_energy,
            'most_stable_energy': most_stable['energy'],
            'most_stable_step': most_stable['step'],
            'energy_improvement': initial_energy - most_stable['energy'],
            'md_parameters': {
                'temperature': temperature,
                'steps': steps,
                'sample_interval': sample_interval
            },
            'energy_analysis': energy_analysis,
            'total_samples': len(sampled_structures),
            'has_surface': has_metal
        }
        
        import json
        with open(f"{output_dir}/md_stability_analysis.json", 'w') as f:
            json.dump(analysis_report, f, indent=2, default=str)
        
        # 保存能量数据
        energy_data = {
            'steps': [s['step'] for s in sampled_structures],
            'energies': [s['energy'] for s in sampled_structures]
        }
        
        with open(f"{output_dir}/energy_data.json", 'w') as f:
            json.dump(energy_data, f, indent=2)
        
        print(f"\n✅ MD stability search complete!")
        print(f"   📁 Results saved in: {output_dir}/")
        print(f"   🎯 Most stable energy: {most_stable['energy']:.4f} eV")
        print(f"   📈 Energy improvement: {initial_energy - most_stable['energy']:.4f} eV")
        print(f"   📊 Total samples: {len(sampled_structures)}")
        
        return {
            'most_stable': most_stable,
            'analysis': analysis_report,
            'sampled_structures': sampled_structures
        }
        
    except Exception as e:
        print(f"❌ Error during MD stability search: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="MD Stability Search Tool")
    parser.add_argument("structure_file", help="Structure file path (CIF, XYZ, etc.)")
    parser.add_argument("--temperature", type=float, default=300, 
                       help="MD simulation temperature in K (default: 300)")
    parser.add_argument("--steps", type=int, default=100000,
                       help="Number of MD steps (default: 100000)")
    parser.add_argument("--sample_interval", type=int, default=100,
                       help="Sampling interval (default: 100)")
    parser.add_argument("--output_dir", type=str, default="md_stability_results",
                       help="Output directory name (default: md_stability_results)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.structure_file):
        print(f"❌ Structure file not found: {args.structure_file}")
        return
    
    print("🧪 MD Stability Searcher")
    print("=" * 50)
    print(f"📖 Structure file: {args.structure_file}")
    print(f"🌡️  Temperature: {args.temperature} K")
    print(f"🔄 Steps: {args.steps}")
    print(f"📊 Sample interval: {args.sample_interval}")
    print(f"📁 Output directory: {args.output_dir}")
    print()
    
    # 运行MD稳定性搜索
    result = run_md_stability_search(
        args.structure_file, 
        temperature=args.temperature,
        steps=args.steps,
        sample_interval=args.sample_interval
    )
    
    if result:
        print("\n🎉 MD stability search completed successfully!")
    else:
        print("\n❌ MD stability search failed!")


if __name__ == "__main__":
    main()
