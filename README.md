# Retro Synthesis

基于 **AiZynthFinder** 的黄酮苷类化合物逆合成分析平台，专注于黄酮类糖苷（如橙皮苷）的逆合成路线规划。

## 项目特点

- 从文献笔记中程序化提取 508 条黄酮逆反应模板（16 类反应）
- 自定义 `UniformTemplateExpansion` 扩展策略，无需训练模型即可使用小型领域模板库
- 支持 USPTO 通用模板 + 黄酮专用模板的联合搜索
- 包含完整的模板提取、生成、转换、验证流水线

## 目录结构

```
retro_synthesis/
├── src/                    # 源代码 (提取/生成/转换/验证流水线)
├── data/                   # 预训练模型 & 股票数据库 (需单独下载)
├── templates/              # 自定义黄酮反应模板库
├── outputs/                # 搜索结果
├── files/                  # 学习笔记 & 改进计划
├── docs/                   # 文档 (SOP, 教程)
├── notebooks/              # Jupyter notebooks
├── config/                 # 配置文件
├── scripts/                # 独立脚本
└── structure.md            # 详细项目结构说明
```

## 环境配置

### 依赖

- Python 3.13 (模板提取/验证)
- Python 3.10-3.12 + AiZynthFinder (逆合成搜索，需 conda 环境 `retro`)
- RDKit, pandas, matplotlib, seaborn, networkx, tables

### 安装

```bash
# 创建虚拟环境
python3 -m venv venv
source venv/bin/activate

# 安装依赖
pip install -r requirements.txt

# AiZynthFinder 需要单独的 conda 环境
conda create -n retro python=3.11
conda activate retro
pip install aizynthfinder
```

### 数据文件

`data/` 目录下的大型文件（模型、股票数据库）未纳入版本控制，需单独下载：

```bash
# USPTO 模型和模板
# ZINC 股票数据库
# 详见 data/README.md (如有)
```

## 快速开始

### 模板提取

```bash
# 从文献笔记提取黄酮反应模板
python src/build_literature_templates.py

# 验证模板
python src/validate_templates_on_hesperidin.py
```

### 逆合成搜索

```bash
# 使用联合配置运行 AiZynthFinder
aizynthfinder --config config/flavonoid.yml --smiles "COc1ccc(...)"

# 可视化结果
python scripts/processing.py
```

## 文档

- [项目结构](structure.md) — 详细的目录和文件说明
- [AiZynthFinder SOP](docs/AIZYNTHFINDER_SOP.md) — 标准操作流程
- [改进计划](files/retrosynthesis_improvement_plan.md) — 12 周优化路线图
- [优化指南](files/基于%20AiZynthFinder%20的橙皮苷类分子逆合成平台优化指南.md)
