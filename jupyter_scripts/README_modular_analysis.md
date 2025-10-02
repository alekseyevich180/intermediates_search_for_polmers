# 模块化分析工具使用说明

## 📋 概述

现在search.py的所有功能已经拆分成独立的脚本，每个都可以单独运行，都接受结构文件作为输入。这样的设计让您可以：

1. **灵活选择** - 只运行需要的分析模块
2. **独立调试** - 每个模块可以单独测试和调试
3. **并行运行** - 可以同时运行多个模块
4. **易于扩展** - 添加新的分析模块很简单

## 🔧 可用的分析模块

### 1. `surface_extractor.py` - 表面分子提取器
**功能**: 从表面+分子结构中提取分子部分
**输入**: 表面+分子结构文件
**输出**: 分离的表面和分子结构

```bash
# 基本用法
python surface_extractor.py your_surface_structure.cif

# 自定义参数
python surface_extractor.py your_file.cif --surface_threshold 0.3 --molecule_threshold 0.4
```

### 2. `functional_group_analyzer.py` - 官能团分析器
**功能**: 识别分子中的官能团和周围的C原子
**输入**: 结构文件 (支持表面+分子或纯分子)
**输出**: 官能团信息、C原子位置、键网络

```bash
# 基本用法
python functional_group_analyzer.py your_structure.cif

# 自定义搜索距离
python functional_group_analyzer.py your_file.cif --max_distance 3.5
```

### 3. `md_stability_searcher.py` - MD稳定性搜索器
**功能**: 通过MD模拟寻找最稳定的结构
**输入**: 结构文件 (支持表面+分子或纯分子)
**输出**: 最稳定结构、采样轨迹、能量分析

```bash
# 基本用法
python md_stability_searcher.py your_structure.cif

# 自定义MD参数
python md_stability_searcher.py your_file.cif --temperature 400 --steps 50000 --sample_interval 50
```

### 4. `carbon_optimizer.py` - C原子优化器
**功能**: 从近到远逐步优化C原子位置
**输入**: 结构文件 (支持表面+分子或纯分子)
**输出**: 优化后的结构、优化历史、能量变化

```bash
# 基本用法
python carbon_optimizer.py your_structure.cif

# 自定义优化参数
python carbon_optimizer.py your_file.cif --opt_trials 150
```

### 5. `reaction_detector.py` - 反应检测器
**功能**: 检测和分析氧化反应中的中间状态
**输入**: 结构文件 (支持表面+分子或纯分子)
**输出**: 反应中间体、键变化分析、反应路径

```bash
# 基本用法
python reaction_detector.py your_structure.cif
```

### 6. `analysis_master.py` - 主控制脚本
**功能**: 统一控制所有分析模块，提供完整的分析流程
**输入**: 结构文件 (支持表面+分子或纯分子)
**输出**: 完整的分析报告和结果文件

```bash
# 运行所有模块
python analysis_master.py your_structure.cif

# 运行特定模块
python analysis_master.py your_file.cif --modules functional_group_analyzer carbon_optimizer

# 创建示例结构
python analysis_master.py --create_sample
```

## 🚀 使用流程

### 方法1: 使用主控制脚本 (推荐)

```bash
# 1. 创建示例结构 (可选)
python analysis_master.py --create_sample

# 2. 运行完整分析
python analysis_master.py sample_surface_with_2_nonanone.cif

# 3. 查看结果
ls complete_analysis_results/
```

### 方法2: 逐步运行各个模块

```bash
# 1. 提取分子部分
python surface_extractor.py your_surface_structure.cif

# 2. 分析官能团
python functional_group_analyzer.py your_structure.cif

# 3. MD稳定性搜索
python md_stability_searcher.py your_structure.cif

# 4. C原子优化
python carbon_optimizer.py your_structure.cif

# 5. 反应检测
python reaction_detector.py your_structure.cif
```

### 方法3: 只运行需要的模块

```bash
# 只进行官能团分析和C原子优化
python functional_group_analyzer.py your_file.cif
python carbon_optimizer.py your_file.cif
```

## 📊 输出文件结构

### 每个模块的输出
```
module_name_results/
├── analysis_report.json          # 分析报告
├── extracted_structures.cif      # 提取的结构 (如果有)
├── optimization_summary.json     # 优化总结 (如果适用)
└── other_output_files...         # 其他输出文件
```

### 完整分析的输出
```
complete_analysis_results/
├── analysis_summary.json                    # 总结报告
├── surface_extractor_results/               # 表面提取结果
├── functional_group_analyzer_results/       # 官能团分析结果
├── md_stability_searcher_results/          # MD稳定性搜索结果
├── carbon_optimizer_results/               # C原子优化结果
└── reaction_detector_results/              # 反应检测结果
```

## 🎯 参数说明

### 通用参数
- `structure_file`: 结构文件路径 (必需)
- `--output_dir`: 输出目录名称 (可选)

### 模块特定参数

#### surface_extractor.py
- `--surface_threshold`: 表面层阈值 (默认: 0.3)
- `--molecule_threshold`: 分子层阈值 (默认: 0.4)

#### functional_group_analyzer.py
- `--max_distance`: 搜索C原子的最大距离 (默认: 3.0 Å)

#### md_stability_searcher.py
- `--temperature`: MD模拟温度 (默认: 300K)
- `--steps`: MD模拟步数 (默认: 100000)
- `--sample_interval`: 采样间隔 (默认: 100)

#### carbon_optimizer.py
- `--opt_trials`: 每个C原子的优化试验次数 (默认: 100)

#### analysis_master.py
- `--modules`: 指定要运行的模块
- `--create_sample`: 创建示例结构文件

## 💡 使用建议

### 1. 首次使用
```bash
# 创建示例结构并运行完整分析
python analysis_master.py --create_sample
python analysis_master.py sample_surface_with_2_nonanone.cif
```

### 2. 处理自己的结构
```bash
# 对于表面+分子结构
python analysis_master.py your_surface_molecule.cif

# 对于纯分子结构
python functional_group_analyzer.py your_molecule.cif
python carbon_optimizer.py your_molecule.cif
```

### 3. 调试特定模块
```bash
# 只运行有问题的模块
python problematic_module.py your_structure.cif --verbose
```

### 4. 批量处理
```bash
# 处理多个结构文件
for file in *.cif; do
    python analysis_master.py "$file" --output_dir "results_${file%.cif}"
done
```

## 🔧 高级用法

### 1. 自定义分析流程
```bash
# 创建自定义分析脚本
cat > custom_analysis.sh << 'EOF'
#!/bin/bash
structure_file=$1

# 1. 提取分子
python surface_extractor.py "$structure_file"

# 2. 分析官能团
python functional_group_analyzer.py "$structure_file"

# 3. 只优化前3个C原子
python carbon_optimizer.py "$structure_file" --opt_trials 50
EOF

chmod +x custom_analysis.sh
./custom_analysis.sh your_structure.cif
```

### 2. 并行运行模块
```bash
# 同时运行多个模块
python functional_group_analyzer.py your_file.cif &
python md_stability_searcher.py your_file.cif &
python reaction_detector.py your_file.cif &
wait
```

### 3. 结果后处理
```bash
# 收集所有结果
python -c "
import json
import glob

results = {}
for result_file in glob.glob('*_results/*.json'):
    with open(result_file) as f:
        results[result_file] = json.load(f)

with open('combined_results.json', 'w') as f:
    json.dump(results, f, indent=2)
"
```

## ❓ 常见问题

### Q1: 模块运行失败
**解决方案**:
- 检查结构文件格式是否正确
- 确保所有依赖包已安装
- 查看错误输出信息
- 尝试减少参数值 (如步数、试验次数)

### Q2: 表面识别不准确
**解决方案**:
- 调整surface_threshold和molecule_threshold参数
- 检查结构的z坐标分布
- 手动检查提取结果

### Q3: 优化时间过长
**解决方案**:
- 减少opt_trials参数
- 减少MD步数
- 使用更小的结构进行测试

### Q4: 内存不足
**解决方案**:
- 减少MD步数
- 减少采样频率
- 分批处理大结构

## 🎉 总结

现在您有了完整的模块化分析工具集！

**主要优势**:
✅ 每个模块都可以独立运行
✅ 支持表面+分子结构自动处理
✅ 灵活的参数配置
✅ 完整的分析报告
✅ 易于扩展和维护

**推荐使用流程**:
1. 使用`analysis_master.py`运行完整分析
2. 根据结果选择需要的模块进行深入分析
3. 使用模块特定参数进行优化调整

开始使用吧！🚀
