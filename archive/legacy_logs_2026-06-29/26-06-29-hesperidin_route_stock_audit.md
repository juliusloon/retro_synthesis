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

---

## 后续操作执行记录（2026-06-29 晚间）

### 1. Strict Buyable Stock 构建

- 脚本: `scripts/build_strict_stock.py`
- 过滤规则: MW > 350 或 重原子数 > 25 的条目从 stock 中排除
- 结果: 297 → 261 条目（排除 36 条复杂黄酮糖苷/天然产物）
- 输出: `templates/strict_buyable_stock_inchikeys.txt`
- 审计报告: `logs/26-06-29-strict_stock_audit.md`

### 2. 模板原子比质量门槛

- 修改: `scripts/custom_expansion.py` 新增 `min_atom_ratio` 参数
- 实现方式: 计算 retro template 前体侧/目标侧的原子数比值，过滤比值过低的"塌缩"模板
- 统计（508 模板）:
  - 0.8-1.2 (balanced): 372
  - 0.5-0.8 (mild collapse): 101
  - 0.3-0.5 (moderate collapse): 4
  - < 0.3 (severe collapse): 1
- 配置: `config/flavonoid.yml` 中三个 custom expansion 均设置 `min_atom_ratio: 0.5`
- 0.5 截止仅过滤 5 个严重塌缩模板，保留化学合理的轻度塌缩模板

### 3. 消融实验结果

5 组实验已完成，结果保存在 `outputs/ablation/`。

| 实验 | 总路线 | 已解决 | 最短步 | 平均步 | 自定义步 | 耗时(s) |
|---|---:|---:|---:|---:|---:|---:|
| baseline_zinc | 10 | 10 | 2 | 5.0 | 0 | 93.5 |
| baseline_strict | 6 | 6 | 16 | 18.3 | 0 | 98.0 |
| flavonoid_zinc | 9 | 9 | 2 | 2.0 | 7 | 14.7 |
| flavonoid_strict | 6 | 6 | 2 | 2.0 | 6 | 66.5 |
| custom_only_strict | 6 | 6 | 2 | 2.0 | 6 | 1.3 |

#### 关键发现

1. **Stock 影响巨大**: Baseline + ZINC 找到 2 步路线，但 Baseline + Strict stock 需要 16-20 步。ZINC 中的复杂黄酮糖苷让 baseline 看起来很容易解决，但这些"起始原料"实际不可购买。

2. **Custom 模板是关键**: 在 strict stock 下，Baseline 需要 16-20 步，而 Custom 模板仅需 2 步。Custom 模板将路线长度缩短了 8-10 倍。

3. **Custom 模板可独立工作**: `custom_only_strict`（无 USPTO）也能找到 6 条 2 步路线，说明黄酮模板自身覆盖了关键逆合成断键。

4. **flavonoid_structural 贡献最大**: 所有 custom 路线均使用 `flavonoid_structural` policy，`flavonoid_biosynthesis` 贡献了 2 条路线。

5. **原子比门槛有效**: 设置 `min_atom_ratio: 0.5` 后，严重塌缩模板被过滤，但保留了化学合理的模板。

6. **无糖苷键断裂**: 所有路线均未检测到糖苷键断裂步骤，这与橙皮苷的结构有关——当前模板库中 `glycosylation` 类模板（6 个）可能需要更多针对橙皮苷糖苷键的逆合成模板。

### 4. 文件变更清单

| 操作 | 文件 | 类型 |
|---|---|---|
| strict stock | `scripts/build_strict_stock.py` | 新建 |
| strict stock | `templates/strict_buyable_stock_inchikeys.txt` | 新建 |
| strict stock | `logs/26-06-29-strict_stock_audit.md` | 新建 |
| 原子比门槛 | `scripts/custom_expansion.py` | 修改 |
| 原子比门槛 | `config/flavonoid.yml` | 修改 |
| 消融实验 | `config/ablation_baseline_zinc.yml` | 新建 |
| 消融实验 | `config/ablation_baseline_strict.yml` | 新建 |
| 消融实验 | `config/ablation_flavonoid_zinc.yml` | 新建 |
| 消融实验 | `config/ablation_flavonoid_strict.yml` | 新建 |
| 消融实验 | `config/ablation_custom_only_strict.yml` | 新建 |
| 消融实验 | `scripts/run_ablation_experiments.py` | 新建 |
| 消融实验 | `scripts/compare_ablation_results.py` | 新建 |
| 消融实验 | `outputs/ablation/*.json` | 新建 |
| 消融实验 | `outputs/ablation/ablation_report.md` | 新建 |

---

## 2026-06-30 修正：solved 统计与 atom-map gate

### 发现的 bug

1. **`count_reactions` 多算一步**: 只要节点有 children 就加 1，把根分子也算成 reaction 步骤。
2. **`is_solved` 未检查 `in_stock`**: 对没有 children 的分子直接返回 `True`，导致 `baseline_strict` 被错误报告为 6 条 solved（实际 `is_solved=False`）。
3. **缺少 map retention 审计**: `min_atom_ratio` 只比较模板 SMARTS 两侧原子数，不能检测应用到橙皮苷后的实际 map retention。6 条 custom 路线全部是 `0/43` map retention（scaffold collapse）。

### 修正后的消融结果

| 实验 | AiZynth solved | map-valid solved | 有效 solved | 最短步 |
|---|---:|---:|---:|---:|
| baseline_zinc | 10 | 9 | 9 | 1 |
| baseline_strict | 0 | 0 | 0 | - |
| flavonoid_zinc | 9 | 2 | 2 | 1 |
| flavonoid_strict | 6 | 0 | 0 | - |
| custom_only_strict | 6 | 0 | 0 | - |

### 修正后的结论

1. **Strict stock 暴露了 baseline 的真实困难度**: ZINC stock 下 baseline 有 9 条有效路线，strict stock 下降为 0。说明 baseline 路线严重依赖 ZINC 中的复杂天然产物作为"可购买"前体。

2. **Custom 模板目前无有效贡献**: `custom_only_strict` 和 `flavonoid_strict` 的 AiZynth solved 是 6，但 map-valid solved 是 0。所有 custom 路线都是 scaffold collapse 型伪路线（0/43 map retention），糖、B 环、取代模式被整体丢弃。

3. **Custom 模板在 ZINC stock 下反而降低了有效路线数**: `flavonoid_zinc` 的有效 solved 从 baseline 的 9 降为 2，说明 custom 模板引入了大量噪声路线。

4. **零 retention 模板清单**: template 28、158、198、199、200、227 在橙皮苷上全部 0% retention，应禁用或重写。

### 已执行的修正

| 修正 | 文件 | 内容 |
|---|---|---|
| count_reactions bug | `scripts/compare_ablation_results.py` | 只统计 `type=='reaction'` 节点 |
| is_solved bug | `scripts/compare_ablation_results.py` | 优先使用顶层 `is_solved` 字段；递归检查要求 `in_stock==True` |
| map retention 审计 | `scripts/compare_ablation_results.py` | 新增 `compute_route_map_retention`，报告三列：AiZynth / map-valid / 有效 |
| map retention 过滤 | `scripts/custom_expansion.py` | 新增 `min_retained_map_ratio` 参数，post-filter 应用反应后的 map retention |
| strict stock 分离 unverified | `scripts/build_strict_stock.py` | 缺少分子信息的 InChIKey 不再"保守保留"，单独输出到 `unverified_stock_inchikeys.txt` |

### 正确的当前状态表述

> Strict stock 暴露了 baseline 在真实库存下无法解决（0 条有效路线）；custom policy 能快速生成 AiZynth solved 候选路线，但这些路线全部未通过 atom-map 守恒审计（0% map retention），属于 scaffold collapse 型伪路线。当前有效 custom solved route 应按 0 处理。下一步需：(1) 禁用零 retention 模板，(2) 重写 O-糖苷键断裂等真实断键模板，(3) 在 `min_retained_map_ratio≥0.8` 门槛下重新消融。

