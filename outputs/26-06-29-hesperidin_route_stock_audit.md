# 橙皮苷逆合成路线 stock 与模板可信度审计

生成日期: 2026-06-29

审计对象:

- `outputs/26-06-29-hesperidin_baseline_uspto.json`
- `outputs/26-06-29-hesperidin_flavonoid_combined.json`
- `config/baseline.yml`
- `config/flavonoid.yml`
- `templates/*_stock_inchikeys_verified.txt`

## 结论摘要

当前两组结果的短路线主要不是搜索深度不足，而是终止条件过宽和部分自定义模板不守恒共同造成的。

1. Baseline 的路线确实全部命中了所选 `zinc` stock，但终止前体里包含大分子黄酮糖苷/甲基化类似物。这在 AiZynthFinder 的 stock 判定上成立，但不一定符合“可作为实际起始原料购买”的实验化学含义。
2. Flavonoid combined 新增的 1 步路线中，路线 3-10 的自定义模板存在严重原子映射问题: 目标分子 43 个映射原子有 0 个在前体侧保留。这些路线应视为模板伪阳性，不应计入优化效果。
3. 三个 verified stock 文件完全相同，均为 176 行: `flavonoid`、`flavonoid_structural`、`flavonoid_biosynthesis` 目前没有实际区分库存集合。
4. Combined 新增路线的多个前体虽然在 `zinc` 中，但没有命中 verified flavonoid stock。因此“加入我的库”的表观优化，混合了 custom policy 命中、ZINC stock 终止和自定义 stock 命中三种因素，需要拆开做消融实验。

## Baseline 路线审计

Baseline 共 9 条路线，所有叶子前体均命中 `zinc_stock_17_04_20.hdf5`。

| 代表前体 | InChIKey | 命中来源 | 备注 |
|---|---|---|---|
| `COc1ccc(... )cc1OC` | GUMSHIGGVOJLBP-FZIAYWRKSA-N | zinc | 橙皮苷甲基化类似物，导致 1 步路线 |
| `C[C@@H]1O[C@@H](OC[C@H]2...` | OMQADRGFMLGFJF-GDRLQRQZSA-N | zinc | 去甲基橙皮苷/相关糖苷大分子，导致 1 步去甲基路线 |
| `COS(=O)(=O)OC` | VAYGXNSJCAHWJZ-UHFFFAOYSA-N | zinc | 小分子试剂，合理 |
| `C=CCBr` | BHELZAPQIKSEDF-UHFFFAOYSA-N | zinc + verified flavonoid stocks | 小分子试剂，合理 |
| `CC(=O)OC(C)=O` | WFDIJRYMOXRFFG-UHFFFAOYSA-N | zinc + verified flavonoid stocks | 小分子试剂，合理 |

Baseline 的 USPTO 反应映射是守恒的: 每一步目标侧映射原子均在前体侧保留，例如 43/43、46/46、52/52。短路线问题主要来自 `zinc` stock 把复杂类似物视为可购买前体。

## Combined 路线审计

Combined 共 11 条路线，其中 1、2、11 为 USPTO 路线，映射守恒；3-10 为自定义 structural/biosynthesis 路线，映射不守恒。

| 路线 | policy | template_code | classification | 映射保留 | 判断 |
|---|---|---:|---|---:|---|
| 1 | uspto | 9071 | Unrecognized | 43/43 | 可作为 stock 审计对象 |
| 2 | uspto | 19597 | Unrecognized | 43/43 | 可作为 stock 审计对象 |
| 3 | flavonoid_structural | 21 | flavonoid_synthesis | 0/43 | 模板伪阳性 |
| 4 | flavonoid_structural | 29 | chalcone_synthesis | 0/43 | 模板伪阳性 |
| 5 | flavonoid_structural | 162 | flavonoid_synthesis | 0/43 | 模板伪阳性 |
| 6 | flavonoid_structural | 202 | chalcone_synthesis | 0/43 | 模板伪阳性 |
| 7 | flavonoid_structural | 203 | chalcone_synthesis | 0/43 | 模板伪阳性 |
| 8 | flavonoid_structural | 204 | flavanone_synthesis | 0/43 | 模板伪阳性 |
| 9 | flavonoid_biosynthesis | 11 | biosynthesis | 0/43 | 模板伪阳性 |
| 10 | flavonoid_biosynthesis | 14 | biosynthesis | 0/43 | 模板伪阳性 |
| 11 | uspto | 28518 | Unrecognized | 43/43 | 可作为 stock 审计对象 |

自定义路线 3-10 的典型症状是: 模板只匹配目标分子的局部芳环或黄酮片段，然后把糖、B 环或其他未映射结构直接丢弃，产物变成一个小的苯乙酮/查尔酮前体。这不是合理逆合成断键。

## Stock 文件状态

当前 verified stock:

| 文件 | 行数 | 关系 |
|---|---:|---|
| `templates/flavonoid_stock_inchikeys_verified.txt` | 176 | 基准 verified 集 |
| `templates/flavonoid_structural_diversity/flavonoid_structural_stock_inchikeys_verified.txt` | 176 | 与基准文件完全相同 |
| `templates/flavonoid_biosynthesis/flavonoid_biosynthesis_stock_inchikeys_verified.txt` | 176 | 与基准文件完全相同 |

原始文献库存候选:

| 文件 | 行数 | 用途建议 |
|---|---:|---|
| `templates/literature_curated/flavonoid_literature_stock_inchikeys.txt` | 297 | 不应直接作为高置信 buyable stock |
| `templates/literature_curated/flavonoid_literature_all_reactants_stock_inchikeys.txt` | 539 | 更适合做模板 reactant 审计，不适合作为默认 stock |

## 后续建议

### 1. 先建立严格可比的消融实验

建议至少跑以下配置矩阵:

| 实验 | expansion | stock | 目的 |
|---|---|---|---|
| baseline_zinc | USPTO | ZINC | 复现默认表现 |
| baseline_strict | USPTO | strict buyable | 测 baseline 在真实库存下的困难度 |
| flavonoid_policy_zinc | USPTO + custom | ZINC | 测 custom policy 是否带来候选步骤 |
| flavonoid_policy_strict | USPTO + custom | strict buyable | 测最终可解决路线是否真实改善 |
| custom_only_strict | custom only | strict buyable | 测黄酮模板自身贡献 |

其中 `strict buyable` 不应包含复杂黄酮糖苷、目标近邻天然产物、只在 ZINC 中出现但不确定可买的高级中间体。

### 2. 给模板加硬性质量门槛

建议在 custom expansion 前或模板构建阶段加入:

- 反应映射覆盖率: target-side mapped atoms retained in precursor side >= 0.8，理想情况下接近 1.0。
- 禁止 0 原子保留模板进入 search。
- 禁止只匹配局部子结构却删除目标其余原子的 structural scaffold 模板。
- 对 `biosynthesis` 模板单独处理: 生物合成关系可作为启发式相似性或 pathway annotation，不应直接当作实验逆合成断键。
- 尽量启用 `rdchiral` 或至少用反应产物原子守恒审计做后验过滤。

### 3. 重新定义“优化”的评价指标

不要只看 solved route 数和最短步数。建议同时报告:

- 有效路线数: 所有步骤均通过映射守恒审计。
- strict-stock solved route 数。
- 每条路线的 stock 来源: `zinc` / `verified_buyable` / `literature_candidate`。
- 自定义模板贡献: custom step 数、custom step 是否守恒、是否有文献来源。
- 路线化学合理性: 是否先断糖苷键、是否能到查尔酮/黄烷酮/苷元等合理黄酮中间体。

### 4. 对黄酮类目标的路线策略

橙皮苷这类黄酮糖苷，最可信的 custom 模板优先级应是:

1. O-糖苷键形成/断裂，生成苷元 + 合理糖基供体。
2. 黄烷酮/查尔酮互变或环化路线。
3. 酚羟基保护/脱保护、甲基化/脱甲基化。
4. 黄酮/黄酮醇/异黄酮骨架互转，但必须保持原子映射和取代模式。

低优先级或暂时禁用:

- 只描述天然产物结构多样性的 scaffold collapse 模板。
- 生物合成相邻物种模板，除非能转化成具体实验反应。
- 任何会把橙皮苷一跳缩成小苯乙酮/查尔酮且无原子保留的模板。

## 当前结果可如何表述

比较稳妥的表述是:

> 已成功将文献整理迁移为 AiZynthFinder 可加载的 templates 和 stock，并能改变橙皮苷搜索返回的候选路线分布。但当前 observed improvement 主要表现为更多 1 步 solved routes；经审计，其中 baseline 短路线来自 ZINC 对复杂类似物的 stock 判定，combined 的多数新增 custom routes 存在原子映射不守恒。因此，下一阶段应先完成 strict stock 和 atom-conserving template gate，再评估黄酮专用模板是否带来真实路线优化。

