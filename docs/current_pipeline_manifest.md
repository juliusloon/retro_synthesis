# 当前逆合成管线 Manifest：黄酮苷模板优化研究架构

最后更新：2026-07-01

## 结论先行

当前 pipeline **还不能算很好地构建了**“使用 AiZynthFinder 的 USPTO/RingBreaker 通用扩展策略 + ZINC 通用库存 + 分层 in-stock 证据 + 人工审核黄酮反应族模板，评估并优化橙皮苷/黄酮苷逆合成”的研究架构。

它已经有正确的骨架：`ZINC baseline`、`strict/trusted`、`virtual_bridge`、人工审核模板、panel 评估、sugar artifact 审计都在场。但它现在更像一个快速迭代后的运行清单，而不是一个能支撑论文式结论的实验架构。最大问题不是“缺脚本”，而是**语义边界和证据层级没有被硬编码成研究问题**。

关键判断：

1. **ZINC 是库存基线，不是反应模板。** 不能把 “USPTO + ZINC + in-stock templates” 混成一个层；必须拆成 expansion policy、stock layer、evidence label、post-hoc audit 四个正交维度。
2. **A1-C2 新框架已在配置和 runner 中出现，但当前输出仍主要是旧命名实验。** Manifest 同时写“已运行”和“待重跑”，会污染结论。
3. **ZINC solved 只能回答大库存可闭合性，不能证明合成可行性。** 需要单独惩罚 target-like natural product stock、保护态 sugar artifact、虚拟 bridge 依赖和芳香黄酮苷叶子。
4. **人工审核模板还没有形成闭环优化。** 现在是“启用/禁用模板 + 事后审计”；研究上更强的结构应是 gap mining -> template proposal -> mapped audit -> positive/negative controls -> ablation -> promotion/demotion。
5. **黄酮苷的真正科学问题应从“是否 solved”升级为“断对了什么键、闭合到了什么证据层”。** 橙皮苷/rutinose、neohesperidose、rhamnoside 的失败和成功机制必须分开报告。

因此，本 manifest 将当前 pipeline 重构为一个更大胆但可审计的架构：**Discovery-Evidence-Mechanism 三轴评估框架**。

## 核心研究命题

主问题：

> 人工审核的黄酮反应族模板，是否在 AiZynthFinder 的 USPTO/RingBreaker 通用策略和 ZINC 库存基线之上，带来可解释、可泛化、可审计的橙皮苷/黄酮苷逆合成改进？

拆成三个可检验命题：

| 命题 | 问题 | 主证据 |
|---|---|---|
| H1 Discovery | 大库存下能否找到闭合路线？ | ZINC solved、route diversity、target-like stock penalty |
| H2 Mechanism | 新模板是否让路线经过合理黄酮断键？ | glycoside cleavage、aglycone/chalcone reach、family knockout |
| H3 Evidence | 路线终点是否具有可购买或可信中间体证据？ | strict/trusted solved、bridge-closed、donor sandbox audit |

任何单一 `is_solved` 都不能回答这三个问题。

## 四层架构

### Layer 1：Expansion Policy 层

只描述“能生成哪些逆合成动作”，不描述库存。

| Policy | 角色 | 进入主结论条件 |
|---|---|---|
| `uspto` | AiZynthFinder 通用模板基线 | 始终作为 A1 基线 |
| `ringbreaker` | 通用开环补充 | 独立消融，不能和 custom 增益混报 |
| `flavonoid_reaction_families` | 人工审核黄酮反应族模板 | 必须通过 map retention、family metadata、negative controls |
| `glycosyl_donor_surrogates` | 糖供体候选模板沙盒 | 默认禁用；只能做 donor identity audit |

当前活跃 custom 模板只有：

- `o_glycoside_cleavage_pyranose`
- `aryl_methyl_ether_cleavage`
- `aryl_O_methylation`

当前应继续禁用：

- chalcone/flavanone 模板：除非重写 SMARTS 并通过实际 target/panel map retention。
- activated donor 模板：除非编码离去基、供体 exact identity 和负对照。
- 任何 scaffold collapse 或 0 map retention 模板。

### Layer 2：Stock/Evidence 层

只描述“叶子节点能否作为终点解释”，不描述模板增益。

| 层级 | 角色 | 可报告标签 | 禁止事项 |
|---|---|---|---|
| `ZINC` | 主 discovery 库存 | `zinc_discovery_solved` | 不能等同可购买合成路线 |
| `ZINC-filtered` | 建议新增：去除 target-like 黄酮苷/复杂天然产物的 ZINC 子集 | `zinc_filtered_solved` | 不能手工挑选有利条目 |
| `strict_buyable` | 高置信可购买小分子 | `strict_solved` | 不接收 virtual bridge 或 donor surrogate |
| `trusted_intermediate` | 文献/家族支持中间体 | `trusted_solved` | 不等同供应商可购买 |
| `virtual_bridge` | 糖层连通性诊断 | `bridge_closed_connectivity` | 不能作为真实 solved |
| `donor_sandbox` | 糖供体 exact identity 审计 | `donor_audit_only` | 默认不进入 MCTS 主结论 |

Bold reconstruction：新增 `ZINC-filtered` 是必要的。当前 ZINC baseline 会把复杂天然产物、目标近邻黄酮苷或保护态糖片段当作终点，这对 discovery 有价值，但对“优化模板”会产生强烈混淆。

### Layer 3：Evaluation 层

所有路线必须先经过共同质量门：

| Gate | 阈值/规则 | 目的 |
|---|---|---|
| atom-map retention | 每步 `>=0.8` | 防 scaffold collapse |
| no target-like leaf | 黄酮苷/目标近邻叶子单独标记 | 防 ZINC 偷懒闭合 |
| protected sugar artifact flag | 保护态 sugar leaf 惩罚 | 防 USPTO 保护基偏置冒充 donor |
| virtual bridge isolation | bridge route 单独分类 | 防诊断库存变成结论库存 |
| exact donor identity | donor sandbox 专用 | 防 neutral sugar/bridge artifact 冒充 activated donor |

推荐路线标签优先级：

1. `strict_solved`
2. `trusted_solved`
3. `zinc_filtered_solved`
4. `zinc_discovery_solved`
5. `bridge_closed_connectivity`
6. `donor_audit_only`
7. `unsolved`

`bridge_closed_connectivity` 和 `donor_audit_only` 是科学诊断，不是合成成功。

### Layer 4：Human Review/Promotion 层

人工审核不是最后写一句“已审核”，而是一个 promotion system。

| 状态 | 含义 | 可进入 |
|---|---|---|
| `proposed` | 从 route gap、文献或机制直觉提出 | worklist |
| `mapped_audit_pass` | SMARTS 可运行且保留 atom map | inactive sandbox |
| `control_pass` | positive target 命中、negative sugar family 不误命中 | ablation candidate |
| `ablation_positive` | A/B 消融带来机制性改进 | active custom policy |
| `evidence_promoted` | stock/donor exact evidence 通过 | trusted/strict candidate |
| `demoted` | artifact、方向敏感、identity 不足 | audit only |

这个 promotion system 应写入模板和库存元数据，避免靠记忆维护。

## 新实验矩阵

### Axis A：Discovery 模板增益，库存固定为 ZINC

| 实验 | Expansion | Stock | 回答 |
|---|---|---|---|
| A1 | USPTO | ZINC | 通用基线 |
| A2 | USPTO + RingBreaker | ZINC | RingBreaker 是否有增益 |
| A3 | USPTO + RingBreaker + Custom | ZINC | 人工黄酮模板是否有增益 |

主比较：`A3 - A2`，不是 `flavonoid_zinc - baseline_zinc` 的旧命名比较。

### Axis B：Evidence 增益，Expansion 固定为 A3

| 实验 | Expansion | Stock | 回答 |
|---|---|---|---|
| B0 | A3 | ZINC-filtered | 去掉 target-like stock 后是否仍可闭合 |
| B1 | A3 | strict | 高置信可购买是否闭合 |
| B2 | A3 | strict + trusted | 文献/中间体证据是否闭合 |
| B3 | A3 | strict + trusted + virtual_bridge | 糖层连通性是否闭合 |

注意：如果实际 AiZynthFinder 需要同时加载 ZINC 才能保持搜索空间，可以运行两个版本：

- `B*_pure`: 只用该证据层，回答真实证据闭合。
- `B*_overlay`: ZINC + 该证据层，回答在 discovery 搜索中哪些叶子证据更硬。

不要把这两种口径混报。

### Axis C：Custom 模板独立性

| 实验 | Expansion | Stock | 回答 |
|---|---|---|---|
| C1 | Custom only | ZINC | custom 是否能独立提出路线 |
| C2 | Custom only | ZINC-filtered | custom 是否依赖 target-like stock |
| C3 | Custom only | strict + trusted + virtual_bridge | custom 是否能独立完成机制诊断 |

Custom-only 的目标不是替代 USPTO，而是验证人工模板是否真能驱动关键断键。

### Axis D：Mechanism family knockout

| 实验 | 改动 | 回答 |
|---|---|---|
| D1 | remove `o_glycoside_cleavage_pyranose` | 糖苷键断裂是否是核心贡献 |
| D2 | remove aryl OMe cleavage/methylation | OMe 步骤是否只是噪声 |
| D3 | add rhamnoside-specific template sandbox | quercitrin 失败是否来自单糖苷断键缺口 |
| D4 | add rewritten chalcone/flavanone sandbox | 骨架断裂是否能带来更深路线 |
| D5 | add donor exact-identity sandbox | activated donor 能否从诊断变成证据 |

### Axis P：Panel 泛化

当前 panel 保留：

- hesperidin：rutinose/flavanone glycoside
- narirutin：rutinose/flavanone glycoside
- rutin：rutinose/flavonol glycoside
- naringin：neohesperidose/flavanone glycoside
- neohesperidin：neohesperidose/flavanone glycoside
- quercitrin：rhamnoside/flavonol glycoside

报告必须按糖家族分层：

| 糖家族 | 主要问题 |
|---|---|
| rutinose | disaccharide bridge 与 exact donor identity |
| neohesperidose | 2-O linkage 是否被误归入 rutinose-like bridge |
| rhamnoside | 单糖苷断键模板缺口 |

## 当前文件状态解读

### 可以信任的部分

| 路径 | 判断 |
|---|---|
| `config/reaction_families.json` | 可作为当前人工模板来源；但需要 promotion metadata 更结构化 |
| `templates/reaction_families/flavonoid_reaction_family_templates.hdf5` | 可作为 active custom policy；只应包含 active templates |
| `templates/reaction_families/flavonoid_reaction_family_templates_audit.csv` | 可作为模板审计表 |
| `templates/stock_layers/*` | 分层库存方向正确；virtual bridge 边界必须继续保持 |
| `config/flavonoid_target_panel.csv` | panel 结构验证方向正确 |
| `outputs/panel_ablation/panel_ablation_report.md` | 可作为当前 panel 现象摘要，但实验名仍是旧三组口径 |

### 必须修正的部分

| 问题 | 风险 | P0 修正 |
|---|---|---|
| manifest 同时声称 A1-C2 活跃且结果快照为旧框架 | 结论不可复现 | 重新跑 A1-C2，并把旧输出归档或标记 legacy |
| `outputs/ablation` 仍是旧命名结果 | gain analysis 与新矩阵不一致 | 生成 `uspto_zinc.json` 等新命名输出 |
| `compare_ablation_results.py` 曾同时含新 gain 和旧变量结论 | 可能运行时报错或写出混合报告 | 已做最小防混修复：识别 `A1-C2`/`legacy` 结果集；重跑后仍应输出 A1-C2-only active report |
| panel 报告仍用 `flavonoid_zinc/flavonoid_virtual_bridge/custom_only_virtual_bridge` | 与 A1/A3/B3 名称不一致 | 重跑或做 migration map |
| 没有 ZINC-filtered 层 | ZINC solved 可能被 target-like stock 高估 | 建立 target-like/natural-product stock filter |
| artifact 惩罚未进入路线排序 | 只能事后解释 | 在 report 中加入 penalized score 和排序 |

## 当前结果只能这样读

旧单靶标结果可以作为历史观察：

- `baseline_zinc` 在橙皮苷上 solved 很强，说明 ZINC 大库存足以闭合，但不能说明路线真实可买。
- `flavonoid_zinc` 比 `baseline_zinc` solved 少，说明 custom policy 可能改变搜索分布并引入竞争，而不是简单增益。
- `flavonoid_virtual_bridge` bridge-closed 强，说明糖层连通性可被诊断闭合，但不等于 strict/trusted solved。

当前 panel 观察可以作为假设生成：

- neohesperidose 靶标在 ZINC 下表现好，可能不是模板瓶颈。
- rutinose 靶标在 ZINC 下部分可闭合，但 bridge identity/donor exactness 仍是核心瓶颈。
- quercitrin/rhamnoside 失败最值得优先攻克，可能代表单糖苷断键模板缺口。

但这些还不能作为新 A1-C2 架构下的最终结论。

## 推荐执行序列

从仓库根目录运行。

### 0. 先冻结旧结果

```bash
mkdir -p archive/legacy_ablation_pre_A1C2_2026-07-01
mv outputs/ablation/baseline_*.json archive/legacy_ablation_pre_A1C2_2026-07-01/
mv outputs/ablation/flavonoid_*.json archive/legacy_ablation_pre_A1C2_2026-07-01/
mv outputs/ablation/custom_only_*.json archive/legacy_ablation_pre_A1C2_2026-07-01/
```

不要删除旧结果；它们对问题发现有价值。

### 1. 重建模板

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_reaction_family_templates.py
```

### 2. 跑新单靶标矩阵

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_ablation_experiments.py --experiments A1,A2,A3,B1,B2,B3,C1,C2
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_ablation_results.py
```

### 3. 跑 panel 核心矩阵

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_panel_ablation_experiments.py --experiments A1,A3,B3
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_panel_ablation_results.py
```

### 4. 跑 gap/artifact/donor 审计

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/generate_route_gap_worklist.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_sugar_bridge_layer.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/audit_sugar_bridge_evidence.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/group_sugar_bridge_cores.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_protected_sugar_artifact_review.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_rutinose_donor_evidence_worklist.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_glycosyl_donor_surrogates.py
```

## Bold reconstruction backlog

### P0：让当前架构可复现

1. 完成 `compare_ablation_results.py` 的结果集识别后，重跑 A1-C2，并让 active report 只承载 A1-C2 主结论；legacy 结果只作历史附录。
2. 重新生成 `outputs/ablation/uspto_zinc.json`、`uspto_rb_zinc.json`、`uspto_custom_zinc.json` 等新框架输出。
3. 把旧 `baseline_zinc/flavonoid_zinc/...` 结果归档为 legacy，不再出现在 active report 中。
4. panel 报告统一使用 A1/A3/B3 名称，保留旧名 migration table。

### P1：加入 ZINC-filtered

建立 `templates/stock_layers/zinc_filtered_exclusion.csv`，过滤：

- 目标分子本身和 stereoisomer/salt。
- 大分子黄酮苷和目标近邻天然产物。
- 已被判定为 protected sugar artifact 的条目。
- 缺供应商证据但高度接近目标的 complex glycoside。

输出两个指标：

- `zinc_discovery_solved`
- `zinc_filtered_solved`

二者差值就是 ZINC 对路线难度的“作弊空间”。

### P2：把 route score 改成 evidence-aware score

建议新增排序分：

```text
route_score =
  aizynth_score
  + mechanism_bonus(glycoside cleavage, aglycone/chalcone reach)
  + evidence_bonus(strict/trusted leaves)
  - artifact_penalty(protected sugar, target-like leaf)
  - bridge_penalty(virtual bridge dependence)
  - depth_penalty(excessive steps)
```

这个分数只用于排序和解释，不替代原始 solved label。

### P3：优先攻克 quercitrin/rhamnoside

quercitrin 是当前最好的 stress test。建议建立 rhamnoside-specific sandbox：

- positive controls：quercitrin、已验证 rhamnoside analogs。
- negative controls：rutinose/neohesperidose disaccharide targets。
- 必须证明它不是把所有 O-glycoside 都粗暴切掉。

### P4：donor exact identity 从 sandbox 升级为证据路线

当前 `C12H22O9` bridge artifact 到 `C24H33ClO15` rutinosyl chloride donor 的连通性映射只能证明 connectivity，不证明 exact stereochemical identity。升级条件：

1. donor 候选结构 exact InChIKey/stereo block 通过。
2. activated leaving group 被模板显式编码。
3. beta/rutinose 与 neohesperidose 负对照均通过。
4. 有供应商、文献实验或可合成前体证据。

满足前不要进入 strict/trusted，只能保持 donor sandbox。

## Active/Legacy 边界

Active：

- `config/ablation_A*.yml`
- `config/ablation_B*.yml`
- `config/ablation_C*.yml`
- `scripts/run_ablation_experiments.py`
- `scripts/run_panel_ablation_experiments.py`
- `config/reaction_families.json`
- `config/flavonoid_target_panel.csv`
- `templates/stock_layers/*`

Legacy：

- `config/ablation_baseline_*.yml`
- `config/ablation_flavonoid_*.yml`
- `config/ablation_custom_only_*.yml`
- `config/flavonoid.yml`
- `config/hesperidin_optimized.yml`
- `templates/flavonoid_templates.*`
- `templates/flavonoid_structural_diversity/`
- `templates/flavonoid_biosynthesis/`
- current old-name JSON outputs under `outputs/ablation/` until rerun

## 最终报告应回答的五句话

最终论文/报告不应只写“solved 多少条”，而要能回答：

1. 在 USPTO/RingBreaker + ZINC discovery 下，橙皮苷和 panel 哪些能闭合？
2. 加入人工黄酮模板后，闭合数、路线机制和路线排序是否真实改善？
3. 去除 target-like ZINC stock 后，改善是否仍存在？
4. strict/trusted/virtual bridge 分别支持哪些证据强度？
5. 失败靶标，尤其 quercitrin/rhamnoside，指向哪类模板或库存缺口？

这个架构下，我们的结论会更锋利：不是“AiZynthFinder 找到路线”，而是“黄酮专用知识在什么条件下、通过什么机制、以什么证据强度改善了逆合成搜索”。
