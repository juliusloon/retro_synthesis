# Retro Synthesis 项目结构

> 基于 AiZynthFinder 的黄酮苷类化合物逆合成分析平台
> 最后更新: 2026-06-28

## 目录总览

```
retro_synthesis/
├── src/                    # 源代码 (Python 脚本)
├── data/                   # 预训练模型 & 股票数据库 (~1.7 GB, .gitignore 排除)
├── templates/              # 自定义黄酮反应模板库 (~8 MB)
│   ├── flavonoid_biosynthesis/
│   ├── flavonoid_structural_diversity/
│   └── literature_curated/
├── outputs/                # 搜索结果 (JSON + 报告)
├── files/                  # 学习笔记 & 改进计划
├── docs/                   # 文档 (SOP, 教程)
├── notebooks/              # Jupyter notebooks
│   └── logs/               # 早期开发 notebooks
├── config/                 # 配置文件
│   ├── flavonoid.yml       # 主配置 (黄酮专用)
│   ├── baseline.yml        # 基线配置 (USPTO only)
│   └── data.yml            # 数据配置
├── scripts/                # 独立脚本
│   ├── custom_expansion.py # 自定义扩展策略
│   └── processing.py       # 结果可视化 & 分析
├── models/                 # [空] 预留模型目录
├── venv/                   # Python 3.13 虚拟环境 (.gitignore 排除)
├── .gitignore              # Git 忽略规则
├── requirements.txt        # Python 依赖
├── README.md               # 项目说明
└── structure.md            # 本文件
```

## 各目录详细说明

### `src/` — 源代码 (15 个脚本)

流水线分为四个阶段:

#### 提取阶段 (Markdown → 原始反应数据)
| 文件 | 行数 | 说明 |
|------|------|------|
| `build_literature_templates.py` | 1841 | 主流水线: 从 Markdown 提取反应、RDKit 验证、通用基团解析、导出 HDF5/CSV/JSON |
| `extract_flavonoid_biosynthesis.py` | 132 | 解析生物合成 Markdown 中的 SMILES |
| `extract_flavonoid_structures.py` | 92 | 解析结构多样性 Markdown 中的 SMILES |

#### 生成阶段 (原始反应 → 模板 JSON)
| 文件 | 行数 | 说明 |
|------|------|------|
| `generate_biosynthesis_templates.py` | 208 | 生物合成反应 → 模板 JSON |
| `generate_structural_templates.py` | 472 | 手工构建结构多样性模板 (12 结构类 + 3 糖苷 + 10 取代基 + 3 异戊烯基环化 + 3 醌/DA) |

#### 转换阶段 (JSON → AiZynthFinder HDF5)
| 文件 | 行数 | 说明 |
|------|------|------|
| `convert_flavonoid_templates.py` | 177 | 黄酮 JSON → HDF5 |
| `convert_biosynthesis_templates.py` | 163 | 生物合成 JSON → HDF5 (含 SMARTS 清洗) |
| `convert_structural_templates.py` | 178 | 结构多样性 JSON → HDF5 |
| `convert_templates_final.py` | 209 | 通用模板转换器 |

#### 验证阶段
| 文件 | 行数 | 说明 |
|------|------|------|
| `validate_templates_on_hesperidin.py` | 242 | 用橙皮苷验证所有模板的子结构匹配 |
| `validate_biosynthesis_templates.py` | 169 | 用 20 个代表性黄酮验证生物合成模板 |
| `validate_structural_templates.py` | 157 | 用 10 个代表性目标验证结构模板 |
| `diagnose_templates.py` | 79 | 诊断模板解析失败 |
| `diagnose_templates_v2.py` | 67 | 改进版诊断 |

#### 可视化
| 文件 | 行数 | 说明 |
|------|------|------|
| `visualize_routes.py` | 49 | 将 AiZynthFinder 输出 JSON 渲染为路线树图 |

### `data/` — 预训练模型 & 数据 (~1.7 GB)

| 文件 | 大小 | 说明 |
|------|------|------|
| `uspto_model.onnx` | 92 MB | USPTO 扩展策略神经网络 |
| `uspto_templates.csv.gz` | 3.3 MB | 42,554 条 USPTO 反应模板 |
| `uspto_ringbreaker_model.onnx` | 15 MB | 开环策略模型 |
| `uspto_ringbreaker_templates.csv.gz` | 374 KB | 开环模板 |
| `uspto_filter_model.onnx` | 17 MB | 反应过滤模型 |
| `full_uspto_03_05_19_rollout_policy.hdf5` | 300 MB | 完整 USPTO rollout 策略 |
| `full_uspto_03_05_19_unique_templates.hdf5` | 45 MB | 完整 USPTO 唯一模板 |
| `zinc_stock_17_04_20.hdf5` | 1.3 GB | ZINC 可购买化合物库 (2020-04) |
| `zinc_stock.hdf5` | symlink | → zinc_stock_17_04_20.hdf5 |
| `config.yml` | 421 B | 基线配置 |

### `templates/` — 自定义黄酮模板库 (~8 MB)

```
templates/
├── flavonoid_templates.hdf5           # 主模板库 (508 条, 16 类反应)
├── flavonoid_templates.csv.gz         # 完整 CSV (含元数据)
├── flavonoid_stock_inchikeys.txt      # 黄酮起始物料 InChIKey
├── flavonoid_biosynthesis/            # 生物合成模板
│   ├── flavonoid_biosynthesis_templates.hdf5
│   └── flavonoid_biosynthesis_stock_inchikeys.txt
├── flavonoid_structural_diversity/    # 结构多样性模板
│   ├── flavonoid_structural_templates.hdf5
│   └── flavonoid_structural_stock_inchikeys.txt
└── literature_curated/                # 文献精炼模板 (24 文件)
    ├── flavonoid_literature_reactions.json
    ├── flavonoid_literature_templates.hdf5
    ├── flavonoid_literature_templates.csv.gz
    ├── flavonoid_literature_templates_aizynthfinder_ready.hdf5
    ├── flavonoid_literature_templates_aizynthfinder_ready.csv.gz
    ├── flavonoid_literature_templates_metadata.csv
    ├── flavonoid_literature_stock_candidates.csv
    ├── flavonoid_literature_stock_inchikeys.txt
    ├── flavonoid_literature_all_reactants_stock_inchikeys.txt
    ├── flavonoid_literature_invalid_reactions.csv
    ├── generic_group_review.{json,csv,md}
    ├── generic_review_worklist.csv
    ├── manual_review_override_{failures.csv,summary.json}
    └── extraction_audit_report.md
```

### `outputs/` — 搜索结果

| 文件 | 说明 |
|------|------|
| `26-06-12-hesperidin_combined_curated.json` | USPTO + 自定义模板联合搜索 (6 条路线, 最短 2 步) |
| `26-06-12-summary_report.md` | 验证 & 搜索总结报告 |
| `26-06-14-hesperidin_uspto_only.json` | 仅 USPTO 搜索 (5 条路线, 4-6 步) |
| `26-06-14-hesperidin_custom_only.json` | 仅自定义模板搜索 (2 条, 未完全求解) |

### `files/` — 学习笔记 & 计划

黄酮化学学习笔记 (中文, 从这些笔记中程序化提取反应模板):
- 黄酮化合物的结构、生物合成、合成、性质、生物活性、波谱学特征
- 黄酮类化合物 & 异黄酮类化合物研究进展
- 黄酮类化合物立体化学
- `retrosynthesis_improvement_plan.md` — 12 周改进路线图
- `基于 AiZynthFinder 的橙皮苷类分子逆合成平台优化指南.md`

### `docs/` — 文档

| 文件 | 说明 |
|------|------|
| `AIZYNTHFINDER_SOP.md` | AiZynthFinder 标准操作流程 (494 行) |
| `IBM_RXN_橙皮苷逆合成教程.md` | IBM RXN API 逆合成教程 |

### `notebooks/` — Jupyter notebooks

| 文件 | 说明 |
|------|------|
| `aizynth_gui.ipynb` | AiZynthFinder 交互式 GUI |
| `aizynth_results_visual.ipynb` | 搜索结果可视化 |
| `rxn4chemistry.ipynb` | IBM RXN for Chemistry API 集成 |
| `RDKit_scratch.ipynb` | RDKit 分子操作实验 |
| `logs/06-07-notebook.ipynb` | 早期开发 notebook (6 月 7 日) |
| `logs/06-08-notebook.ipynb` | 早期开发 notebook (6 月 8 日) |

### `config/` — 配置文件

| 文件 | 说明 |
|------|------|
| `flavonoid.yml` | 主配置: 5 扩展策略 + 4 股票源 + MCTS 参数 |
| `baseline.yml` | 基线配置: 仅 USPTO + ZINC |
| `data.yml` | 数据目录配置 |

### `scripts/` — 独立脚本

| 文件 | 说明 |
|------|------|
| `custom_expansion.py` | `UniformTemplateExpansion` — 自定义扩展策略 (模板等概率 1.0) |
| `processing.py` | 结果可视化: 统计图表 + 分子结构 + 路线树 + Markdown 报告 |

### 根目录其他文件

| 文件 | 说明 |
|------|------|
| `data.xlsx` | 补充数据表 |
| `.gitignore` | Git 忽略规则 |
| `requirements.txt` | Python 依赖 |
| `README.md` | 项目说明 |

## 关键技术参数

| 参数 | 值 |
|------|-----|
| 目标分子 | 橙皮苷 (Hesperidin), MW ~610 |
| 模板库规模 | 508 条黄酮逆反应模板, 16 类反应 |
| 搜索算法 | MCTS (C=1.4, 500 迭代, 600s, max 10 transforms) |
| 股票库 | ZINC 2020 (1.3 GB) + 黄酮专用 InChIKey 列表 |
| 扩展策略 | `UniformTemplateExpansion` (模板等概率, 无需训练模型) |
| Python 环境 | venv (3.13) + conda `retro` (3.10-3.12, 运行 AiZynthFinder) |
