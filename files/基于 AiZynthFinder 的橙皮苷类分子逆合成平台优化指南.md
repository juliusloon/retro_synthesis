
> 本文档系统阐述如何利用 aizynhfinder 构建一个针对**黄酮-O-糖苷类天然产物**（以[[橙皮苷]] Hesperidin 为代表）的高效逆合成平台。
>
> 核心优化方向：**文献模板库构建** + **聚焦键策略** + **专用 Stock 库** + **多目标评分体系**

---

## 目录

- [[#1. 橙皮苷类分子的逆合成挑战]]
- [[#2. aizynthfinder 核心架构回顾]]
- [[#3. Step 1：构建专用文献模板库]]
- [[#4. Step 2：构建专用 Stock 前体库]]
- [[#5. Step 3：聚焦键策略设计]]
- [[#6. Step 4：Expansion Policy 配置]]
- [[#7. Step 5：Filter Policy 配置]]
- [[#8. Step 6：多目标 Scorer 体系]]
- [[#9. Step 7：搜索参数调优]]
- [[#10. 完整配置模板]]
- [[#11. 实施路线图与检查清单]]
- [[#12. 局限性与人工后处理]]


---

## 1. 橙皮苷类分子的逆合成挑战

### 1.1 [[橙皮苷]]的结构特征

[[橙皮苷]]（Hesperidin）属于**黄酮-O-糖苷类**天然产物，结构核心包括：

| 结构单元 | 特征 | 逆合成难点 |
|---------|------|-----------|
| **黄酮骨架**（C6-C3-C6） | 2-苯基色原酮母核 | 苯环取代基的定位选择性 |
| **O-糖苷键** | 橙皮素-7-O-芸香糖（鼠李糖+葡萄糖） | 糖苷键的温和断裂/构建 |
| **糖链修饰** | 芸香糖 = 鼠李糖(α-1,6)葡萄糖 | 1,6-糖苷键的区域选择性 |
| **甲氧基/酚羟基** | C4'-甲氧基、C5/C7-酚羟基 | 脱甲基/保护基策略 |
| **手性中心** | 糖环多个手性中心、C2-C3 饱和 | 立体选择性控制 |

### 1.2 通用模板库的局限性

aizynthfinder 默认使用的 **USPTO 模板库** 对以下反应类型覆盖**严重不足**：

```
USPTO 覆盖度评估
├── 一般有机反应（酰化、缩合、取代）    ████████████████████ 高
├── 芳香亲电取代 / 傅-克反应           ████████████████░░░░ 中高
├── 黄酮骨架环化（查尔酮→黄酮）        ██████████░░░░░░░░░░ 中
├── 糖基化反应（O-糖苷形成）           ██████░░░░░░░░░░░░░░ 低 ⚠️
├── 糖苷键选择性断裂                   ████░░░░░░░░░░░░░░░░ 低 ⚠️
├── 1,6-糖苷键构建                     ███░░░░░░░░░░░░░░░░░ 极低 ⚠️
└── 保护基正交策略（糖羟基分级保护）     ██░░░░░░░░░░░░░░░░░░ 极低 ⚠️
```

**结论**：必须构建**黄酮-糖苷专用模板库**作为补充。

---

## 2. aizynthfinder 核心架构回顾

在深入优化前，回顾四个可控模块：

```
输入：目标 SMILES
    │
    ├──▶ Expansion Policy（扩展策略）
    │       ├── 模板库（SMARTS + 元数据）
    │       └── 神经网络模型（预测模板概率）
    │
    ├──▶ Filter Policy（过滤策略）
    │       ├── 可行性模型（Keras/ONNX）
    │       ├── BondFilter（freeze_bonds）
    │       └── 自定义过滤器
    │
    ├──▶ Stock（原料库）
    │       └── 可购买/可获取的分子集合
    │
    ├──▶ Scorer（评分器）
    │       ├── StateScorer（默认）
    │       ├── BrokenBondsScorer（聚焦键）
    │       └── RouteCostScorer / NumberOfReactionsScorer...
    │
    └──▶ Search Config（搜索配置）
            ├── max_transforms（最大逆推深度）
            ├── break_bonds / freeze_bonds
            └── algorithm_config（MCTS 参数）
```

---

## 3. Step 1：构建专用文献模板库

### 3.1 需要补充的模板类型

针对橙皮苷类分子，模板库应至少包含以下**六大类反应**：

| 类别 | 反应描述 | SMARTS 示例方向 |
|------|---------|----------------|
| **A. 糖苷化反应** | 酚羟基 + 糖基供体 → O-糖苷 | 三氯乙酰亚胺酯糖基化、Koenigs-Knorr |
| **B. 糖苷键断裂** | O-糖苷 → 糖基 + 苷元 | 酸水解、乙酰解、酶水解 |
| **C. 黄酮骨架构建** | 查尔酮 → 黄酮 | 查尔酮氧化环化 |
| **D. 脱甲基/甲氧基化** | Ar-OCH₃ → Ar-OH | BBr₃ 脱甲基、酶法 |
| **E. 糖链编辑** | 葡萄糖 → 芸香糖 | 1,6-糖基转移、鼠李糖引入 |
| **F. 保护基操作** | 糖羟基保护/脱保护 | 乙酰基、苯甲酰基、硅醚 |

### 3.2 从文献提取 SMARTS 模板

```python
# 示例：将文献反应转化为 SMARTS 模板
from rdkit import Chem
from rdkit.Chem import AllChem

# 文献反应：橙皮素 + 三氯乙酰亚胺酯葡萄糖 → 橙皮苷-7-O-葡萄糖
reaction_smarts = "[OH:1][c:2]>>[O:1]([C@H:3]1O[C@@H:4](CO)[C@H:5](O)[C@@H:6](O)[C@H:7]1O)[c:2]"

# 验证 SMARTS 是否合法
rxn = AllChem.ReactionFromSmarts(reaction_smarts)
assert rxn is not None, "SMARTS 解析失败"
```

**模板数据格式**：

```csv
index,retro_template,template_set,classification,reaction_hash
0,[#8:1]-[C@@H:2]-[C:3]>>[#8:1]-[c:4].[C@@H:2]-[C:3],flavonoid_glycosylation,A.糖苷化,hash0
1,[#8:1]-[c:2]>>[#8:1]-[C@@H:3]1-[C:4]-[C:5]-[C:6]-[C:7]-1,flavonoid_glycosylation,B.糖苷键断裂,hash1
2,[c:1]1[o:2][c:3][c:4][c:5][c:6]1>>[O:7]=[C:8]-[c:9].[c:1]1[c:10]...[c:2]1,flavonoid_cyclization,C.黄酮环化,hash2
```

### 3.3 模板库的组织建议

```
templates/
├── uspto_general.csv.gz           # 通用反应（原有）
├── flavonoid_core.csv             # 黄酮骨架构建
├── glycosylation.csv              # 糖基化反应
├── glycoside_hydrolysis.csv       # 糖苷键断裂
├── deprotection.csv               # 保护基策略
└── custom_mask.json               # 可选：模板掩码
```

**MultiExpansionStrategy 组合**：

```yaml
expansion:
  uspto:
    - models/uspto_keras_model.hdf5
    - templates/uspto_general.csv.gz
  flavonoid:
    - models/flavonoid_keras_model.hdf5
    - templates/flavonoid_core.csv
  glyco:
    - models/glyco_keras_model.hdf5
    - templates/glycosylation.csv
  multi:
    type: aizynthfinder.context.policy.MultiExpansionStrategy
    expansion_strategies: [uspto, flavonoid, glyco]
    additive_expansion: true      # 累加所有策略的结果
    cutoff_number: 80             # 总体截断 Top 80
```

---

## 4. Step 2：构建专用 Stock 前体库

### 4.1 橙皮苷逆合成需要的 Stock 分类

Stock 决定搜索何时终止。如果缺少关键前体，系统会无限拆下去或给出无意义路线。

| 类别 | 示例分子 | SMILES 示例 | 重要性 |
|------|---------|-------------|--------|
| **糖基供体** | 葡萄糖三氯乙酰亚胺酯 | `C1[C@H](OC(=N)C(Cl)(Cl)Cl)...` | ⭐⭐⭐ |
| **鼠李糖供体** | 鼠李糖三氯乙酰亚胺酯 | `C[C@@H]1O[C@@H](OC(=N)C(Cl)(Cl)Cl)...` | ⭐⭐⭐ |
| **黄酮前体** | 橙皮素 | `COc1ccc(C2CC(=O)c3c(cc(O)cc3O)O2)...` | ⭐⭐⭐ |
| **查尔酮前体** | 2,4,6-三羟基苯乙酮 | `CC(=O)c1c(O)cc(O)cc1O` | ⭐⭐⭐ |
| **芳香醛** | 对甲氧基苯甲醛 | `COc1ccc(C=O)cc1` | ⭐⭐ |
| **糖砌块** | 全乙酰葡萄糖 | `CC(=O)OC1C(OC(C)=O)C(OC(C)=O)...` | ⭐⭐ |
| **保护基试剂** | 苯甲酰氯、TBSCl | `O=C(Cl)c1ccccc1` | ⭐⭐ |
| **通用试剂** | 丙二酸、乙酸酐等 | - | ⭐ |

### 4.2 制作 Stock

```bash
# 方式 1：从 SMILES 列表制作 HDF5 stock
smiles2stock \
  --smiles my_stock.smi \
  --output hesperidin_stock.hdf5

# 方式 2：直接提供 CSV（InChI Key + SMILES）
```

```yaml
stock:
  commercial:
    - stocks/zinc_stock.hdf5
  flavonoid_precursors:
    - stocks/flavonoid_precursors.hdf5
  glycosyl_donors:
    - stocks/glycosyl_donors.hdf5
```

---

## 5. Step 3：聚焦键策略设计

### 5.1 橙皮苷分子中的战略键分析

以橙皮苷（Hesperidin）为例，关键键的优先级：

```
              鼠李糖 (Rha)
                 │
          O──C1 (1,6-糖苷键)  ←── break_bonds (高优先级)
                 │
           葡萄糖 (Glc)
                 │
          O──C7 (酚-O-糖苷键)  ←── break_bonds (最高优先级)
                 │
    ┌────────────┴────────────┐
    │      橙皮素 (Hesperetin)  │
    │    ┌──────────────┐     │
    │    │   C4'-OCH₃   │     │
    │    │   (甲氧基)    │     │
    │    └──────────────┘     │
    │         │               │
    │    C3-C4 键？           │
    └─────────────────────────┘
```

### 5.2 键索引的获取方法

```python
from rdkit import Chem
from aizynthfinder.utils.bonds import find_bonds

smiles = "COc1ccc(C2CC(=O)c3c(cc(OC4OC(CO)C(O)C(O)C4O)c(OC4OC(C)C(O)C(O)C4O)c3O)O2)cc1"
mol = Chem.MolFromSmiles(smiles)

# 查看所有原子索引
for atom in mol.GetAtoms():
    print(f"Atom {atom.GetIdx()}: {atom.GetSymbol()}")

# 找到糖苷键的 atom pair（需要结合结构可视化判断）
# 示例：假设 O(糖苷氧) = atom 20, C(黄酮C7) = atom 15
```

### 5.3 freeze_bonds — 保护骨架

```yaml
search:
  freeze_bonds:
    # 保护黄酮母核的 C-C 键，防止把苯环拆开
    - [5, 6]
    - [6, 7]
    - [7, 8]
    # 保护糖环内部 C-C 键（不想拆糖环本身）
    - [25, 26]
    - [26, 27]
    # 保护芳醚 C-O 键（甲氧基）
    - [2, 3]
```

**底层机制**：`BondFilter` 检查反应前后这些键是否断裂，若断裂则 `RejectionException`。

### 5.4 break_bonds — 优先断裂目标键

```yaml
search:
  break_bonds:
    # 最高优先级：酚-O-糖苷键
    - [15, 20]
    # 次优先级：鼠李糖-葡萄糖 1,6-糖苷键
    - [30, 35]
    # 可选：C4'-甲氧基键（若策略是先去甲基）
    - [3, 4]
  algorithm_config:
    search_rewards: ["state score", "broken bonds"]
    search_rewards_weights: [0.7, 0.3]   # 加权模式
    # 或 search_rewards_weights: []      # Pareto 多目标模式
```

**底层机制**：`BrokenBondsScorer` 计算路径中实际断裂了多少目标键，作为额外奖励。

---

## 6. Step 4：Expansion Policy 配置

### 6.1 模型选择

| 场景 | 推荐模型 | 理由 |
|------|---------|------|
| 快速原型 | 预训练 ONNX（USPTO） | 快、无需 TF |
| 专用优化 | 微调 Keras 模型 | 可融入专用模板库 |
| 多策略混合 | MultiExpansionStrategy | 组合通用 + 专用 |

### 6.2 截断参数

```yaml
expansion:
  uspto:
    type: template-based
    model: models/uspto_model.hdf5
    template: templates/uspto_general.csv.gz
    cutoff_cumulative: 0.995
    cutoff_number: 50
    use_rdchiral: true        # 强烈推荐 RDChiral 处理手性
    rescale_prior: false
    chiral_fingerprints: true  # 橙皮苷有手性糖环，建议开启
```

---

## 7. Step 5：Filter Policy 配置

### 7.1 基础过滤

```yaml
filter:
  uspto_filter:
    type: quick-filter
    model: models/uspto_filter.hdf5
    filter_cutoff: 0.05       # 可行性阈值
```

### 7.2 聚焦键过滤（自动注入）

当配置了 `search.freeze_bonds` 时，`AiZynthFinder._setup_focussed_bonds()` 会自动：

1. 创建 `BondFilter` 实例
2. 加载到 `filter_policy`
3. 激活该过滤器

**无需手动配置**，但需确保 YAML 中的 `freeze_bonds` 键对正确。

### 7.3 自定义过滤（可选）

如果你希望过滤掉某些特定反应类别，可以编写自定义 `FilterStrategy`：

```python
# custom_filters.py
from aizynthfinder.context.policy.filter_strategies import FilterStrategy
from aizynthfinder.utils.exceptions import RejectionException

class GlycosideSelectivityFilter(FilterStrategy):
    """只保留 β-构型糖苷化反应的过滤器"""
    
    def apply(self, reaction):
        if "glycosylation" in reaction.metadata.get("template_code", ""):
            if reaction.metadata.get("stereochemistry") != "beta":
                raise RejectionException("非 β-构型糖苷化")
```

---

## 8. Step 6：多目标 Scorer 体系

### 8.1 推荐的 Scorer 组合

对于橙皮苷类分子，建议采用 **加权多目标**：

```yaml
scorer:
  # 默认 scorer，不需要显式配置
```

```python
# 在 search.algorithm_config 中配置奖励
search:
  algorithm_config:
    search_rewards: ["state score", "broken bonds", "number of reactions"]
    search_rewards_weights: [0.5, 0.3, 0.2]
```

| 目标 | 权重建议 | 说明 |
|------|---------|------|
| `state score` | 0.4-0.6 | 核心目标：stock 覆盖率高、步数合理 |
| `broken bonds` | 0.2-0.3 | 战略键断裂奖励 |
| `number of reactions` | 0.1-0.2 | 偏好短路线 |
| `route cost` | 0.1 | 若 stock 有价格数据 |

### 8.2 后处理 Scorer

```yaml
post_processing:
  route_scorers: ["state score", "number of reactions"]
  scorer_weights: [0.6, 0.4]
  min_routes: 10
  max_routes: 50
  all_routes: false
```

---

## 9. Step 7：搜索参数调优

### 9.1 MCTS 核心参数

```yaml
search:
  max_transforms: 8            # 黄酮-糖苷建议 6-10 步
  iteration_limit: 500         # 默认 100 偏少，建议 300-1000
  time_limit: 300              # 5 分钟，复杂分子可延长
  return_first: false          # 不要找到第一个解就停
  algorithm: mcts
  algorithm_config:
    C: 1.4                     # 探索常数，默认即可
    default_prior: 0.5
    use_prior: true            # 使用神经网络的先验概率初始化
    prune_cycles_in_search: true   # 必须开启，防止糖环回环
    immediate_instantiation: []    # 除非调试，不建议开启
```

### 9.2 参数调优建议

| 问题 | 调整方向 |
|------|---------|
| 路线太短、拆不彻底 | 增大 `max_transforms` |
| 路线太长、步数过多 | 增加 `number of reactions` scorer 权重 |
| 找不到糖苷键断裂 | 提高 `broken bonds` 权重，检查 break_bonds 键对 |
| 反复生成相同分子 | 确认 `prune_cycles_in_search: true` |
| 搜索太慢 | 减小 `iteration_limit`，使用 ONNX 模型 |

---

## 10. 完整配置模板

```yaml
# hesperidin_platform_config.yaml

search:
  algorithm: mcts
  max_transforms: 8
  iteration_limit: 500
  time_limit: 300
  return_first: false
  exclude_target_from_stock: true
  
  # 聚焦键策略 —— 根据实际分子调整 atom indices
  break_bonds:
    - [15, 20]    # 酚-O-糖苷键（橙皮素-葡萄糖）
    - [30, 35]    # 1,6-糖苷键（葡萄糖-鼠李糖）
  
  freeze_bonds:
    # 保护黄酮芳香骨架（需根据实际原子索引调整）
    - [5, 6]
    - [6, 7]
    - [7, 8]
    # 保护甲氧基 C-O 键（可选）
    - [2, 3]
  
  algorithm_config:
    C: 1.4
    default_prior: 0.5
    use_prior: true
    prune_cycles_in_search: true
    immediate_instantiation: ()
    mcts_grouping: "partial"
    search_rewards: ["state score", "broken bonds", "number of reactions"]
    search_rewards_weights: [0.5, 0.3, 0.2]

post_processing:
  min_routes: 10
  max_routes: 50
  all_routes: false
  route_scorers: ["state score", "number of reactions"]
  scorer_weights: [0.6, 0.4]

expansion:
  uspto:
    type: template-based
    model: models/uspto_model.onnx
    template: templates/uspto_general.csv.gz
    cutoff_cumulative: 0.995
    cutoff_number: 50
    use_rdchiral: true
    chiral_fingerprints: true
  
  flavonoid:
    type: template-based
    model: models/flavonoid_model.hdf5
    template: templates/flavonoid_core.csv
    cutoff_cumulative: 0.99
    cutoff_number: 30
    use_rdchiral: true
  
  glyco:
    type: template-based
    model: models/glyco_model.hdf5
    template: templates/glycosylation.csv
    cutoff_cumulative: 0.99
    cutoff_number: 30
    use_rdchiral: true
  
  multi:
    type: aizynthfinder.context.policy.MultiExpansionStrategy
    expansion_strategies: [uspto, flavonoid, glyco]
    expansion_strategy_weights: [0.3, 0.3, 0.4]
    additive_expansion: true
    cutoff_number: 80

filter:
  uspto_filter:
    type: quick-filter
    model: models/uspto_filter.hdf5
    filter_cutoff: 0.05

stock:
  commercial:
    - stocks/zinc_stock.hdf5
  flavonoid:
    - stocks/flavonoid_precursors.hdf5
  glycosyl:
    - stocks/glycosyl_donors.hdf5
```

---

## 11. 实施路线图与检查清单

### Phase 1：数据准备（2-4 周）

- [ ] 收集橙皮苷/黄酮-糖苷类文献（糖基化、脱保护、黄酮环化）
- [ ] 从文献反应提取 SMARTS 模板（建议 50-200 条高质量模板）
- [ ] 整理 Stock 库（糖基供体、黄酮前体、保护基试剂）
- [ ] 训练或微调 Expansion 模型（可选，可用预训练模型 + 专用模板库）

| 文件                                                   | 格式         | 内容                      |
| ---------------------------------------------------- | ---------- | ----------------------- |
| `flavonoid_glycoside_reactions.csv` / `.json`        | CSV + JSON | **文献反应数据库**（36条反应记录）    |
| `flavonoid_glycoside_smarts_templates.csv` / `.json` | CSV + JSON | **SMARTS 模板库**（20条精选模板） |
| `flavonoid_glycoside_stock.csv` / `.json`            | CSV + JSON | **Stock 化合物库**（32个化合物）  |
| `flavonoid_project_stats.json`                       | JSON       | 项目统计摘要                  |

#### 数据库内容统计

文献反应数据库（36条）：
- **糖基化反应**：10 条（Schmidt TCA、硫苷、NPTFA、Amberlyst 15、有机催化、Er(OTf)₃、区域选择性 C-7/C-3、二糖供体）
- **脱保护反应**：6 条（Zemplén 脱乙酰、苄基氢解、硅醚脱保护、甲醚脱保护、选择性脱保护、缩酮脱保护）
- **黄酮环化/转化**：10 条（查尔酮→黄烷酮、查尔酮→黄酮、黄烷酮→黄酮、Claisen-Schmidt、Baker-Venkataraman、查尔酮还原、黄烷酮还原、黄酮醇氧化）
- **保护基反应**：5 条（苄基化、全乙酰化、TCA 形成、硅醚化、甲醚化）
- **目标分子**：5 条（橙皮苷、柚皮苷、芦丁、根皮素、根皮苷）

SMARTS 模板库（20条）：按反应类别分布
- **环化/缩合/氧化**：8 条（查尔酮环化、氧化环化、脱氢、Baker-Venkataraman 等）
- **糖基化**：4 条（TCA 供体、硫苷、NPTFA、二糖供体）
- **脱保护**：4 条（乙酰基、苄基、硅醚、甲醚）
- **保护基**：4 条（苄基化、乙酰化、TCA 形成、甲醚化）
Stock 库（32个化合物）
- **糖基供体**：6 个（葡萄糖/半乳糖/甘露糖/鼠李糖的 TCA 供体和硫苷）
- **黄酮前体/苷元**：10 个（橙皮素、柚皮素、槲皮素、根皮素、芹菜素、染料木素等 + A/B 环砌块）
- **保护基试剂**：5 个（溴化苄、乙酸酐、TBSCl、碘甲烷、三氯乙腈）
- **催化剂/试剂**：8 个（TMSOTf、BF₃·OEt₂、NaOMe、TBAF、NIS、DDQ、Pd/C、BBr₃）
- **溶剂**：3 个（DCM、MeOH、AcOH）

### Phase 2：聚焦键分析（1 周）

- [ ] 用 RDKit/Marvin 标注橙皮苷分子的原子索引
- [ ] 确定 `break_bonds` 键对（糖苷键为核心）
- [ ] 确定 `freeze_bonds` 键对（保护黄酮骨架）
- [ ] 验证键对是否存在于目标分子中（`has_all_focussed_bonds`）

### Phase 3：配置与测试（1-2 周）

- [ ] 编写 `config.yaml`（参考第 10 节模板）
- [ ] 单分子测试（橙皮苷 SMILES）
- [ ] 调参：`max_transforms`、`iteration_limit`、scorer 权重
- [ ] 检查路线合理性：糖苷键是否正确断开？stock 前体是否合理？

### Phase 4：批量与扩展（持续）

- [ ] 批量测试类橙皮苷分子（不同糖基取代、不同羟基模式）
- [ ] 路线聚类分析（`RouteCollection.cluster()`）
- [ ] 根据结果迭代补充模板库

---

## 12. 局限性与人工后处理

即使平台优化完善，以下问题仍需要**化学家的专业判断**：

| 局限 | 说明 | 应对策略 |
|------|------|---------|
| **无汇聚式策略** | aizynthfinder 只做线性逆合成 | 人工判断：糖基供体和黄酮骨架可分别合成，最后偶联 |
| **立体化学弱** | 不会自动考虑 α/β 构型 | 人工检查糖苷键构型，必要时加自定义 Filter |
| **保护基策略缺失** | 系统不会建议"先保护再反应" | 在 Stock 中加入保护基修饰的前体，或在模板中编码保护基步骤 |
| **区域选择性** | 多个酚羟基时不知道糖基化位点 | 在 `break_bonds` 中明确指定目标 O-C 键 |
| **反应条件** | 不输出温度、溶剂、催化剂 | 结合文献模板元数据，或对接反应条件预测工具 |
| **生物合成路径** | 完全基于有机化学反应，不包含酶催化 | 可补充酶催化模板（如糖基转移酶）到模板库 |

---

## 附录：橙皮苷参考 SMILES

```text
# 橙皮苷 (Hesperidin)
COc1ccc(C2CC(=O)c3c(cc(OC4OC(COC5OC(C)C(O)C(O)C5O)C(O)C(O)C4O)c(OC)c3O)O2)cc1

# 橙皮素 (Hesperetin，去糖基)
COc1ccc(C2CC(=O)c3c(cc(O)c(OC)c3O)O2)cc1

# 芸香糖 (Rutinose)
OC1OC(COC2OC(C)C(O)C(O)C2O)C(O)C(O)C1O
```

---

> **最后提醒**：aizynthfinder 的逆合成质量上限 = 模板库质量 × Stock 覆盖度 × 聚焦键准确性。对于橙皮苷这类天然产物，**高质量专用模板库**比调参更重要。
