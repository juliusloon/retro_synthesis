# 当前逆合成管线清单

最后更新：2026-06-30

本清单定义了当前活跃的橙皮苷/黄酮苷逆合成工作流。目的是防止旧版模板实验、虚拟桥接假设和严格库存证据混在一起。

## 当前目标

使用 AiZynthFinder 加上已审核的黄酮反应族模板，评估橙皮苷/黄酮苷逆合成，并把糖桥接诊断逐步升级为可审计的糖供体/糖中间体证据。

当前的科学区分如下：

- **strict/trusted solved**（严格/可信已解决）：有意义的合成证据候选。
- **bridge-closed solved**（桥接已关闭已解决）：使用 `virtual_bridge` 进行连通性验证；用于诊断，而非可购买性证明。
- **ZINC solved**（ZINC 已解决）：有用的基线，通常通过将复杂数据库分子视为库存来关闭。
- **donor surrogate sandbox**（糖供体替代沙盒）：未映射、未启用的候选糖供体模板区；当前仅用于记录未来方向，不进入 MCTS 结论。

## 活跃管线

从仓库根目录运行：

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_reaction_family_templates.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_ablation_experiments.py --experiments 1,2,3,4,5,6,7,8,9
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_ablation_results.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_sugar_bridge_layer.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/audit_sugar_bridge_evidence.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/group_sugar_bridge_cores.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/reconstruct_rutinosyl_chloride_candidates.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/validate_rutinosyl_chloride_bridge_mapping.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_rutinose_donor_evidence_worklist.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_glycosyl_donor_surrogates.py
```

跨黄酮 panel 结构已配置并验证，但 panel 搜索不是当前主报告的一部分。需要时运行：

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_panel_ablation_experiments.py --experiments 8,9
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_panel_ablation_results.py
```

## 活跃脚本

| 路径 | 状态 | 作用 |
|---|---|---|
| `scripts/custom_expansion.py` | 活跃 | 统一的自定义 AiZynthFinder 扩展策略，带活跃模板过滤和映射保留后过滤。 |
| `scripts/build_reaction_family_templates.py` | 活跃 | 从 `config/reaction_families.json` 构建已审核的反应族 HDF5/CSV。只有活跃模板会被写入搜索库。 |
| `scripts/run_ablation_experiments.py` | 活跃 | 运行九个实验的库存/模板消融矩阵。 |
| `scripts/compare_ablation_results.py` | 活跃 | 计算 AiZynth 已解决、映射有效已解决、桥接已关闭已解决、库存层诊断和族使用情况。 |
| `scripts/generate_route_gap_worklist.py` | 活跃 | 在消融运行后生成未解析的叶子工作列表。 |
| `scripts/build_sugar_bridge_layer.py` | 活跃 | 构建间隙衍生的非芳香糖虚拟桥接条目和审核表。 |
| `scripts/audit_sugar_bridge_evidence.py` | 活跃审核 | 对 17 个 sugar bridge 条目进行保护态、证据等级和保守决策分级。 |
| `scripts/group_sugar_bridge_cores.py` | 活跃审核 | 去乙酰化保护态 sugar bridge，并按人工确认的糖芯归一化分组。 |
| `scripts/reconstruct_rutinosyl_chloride_candidates.py` | 活跃审核/沙盒 | 从 true rutinose cyclic seed 重建 hexa-O-acetyl-rutinosyl chloride 两个异头候选结构。 |
| `scripts/validate_rutinosyl_chloride_bridge_mapping.py` | 活跃审核/沙盒 | 验证重建 chloride donor 候选去乙酰/去 Cl 补 H 后与 `C12H22O9` bridge artifact 的连通性关系。 |
| `scripts/build_rutinose_donor_evidence_worklist.py` | 活跃审核 | 建立 rutinose donor 证据工作列表，区分命名文献 donor、可机器验证结构和不可升级候选。 |
| `scripts/fix_disaccharide_linkages.py` | 活跃审核/一次性修正 | 区分 rutinose 与 neohesperidose 的 6-O/2-O 连接语义，并更新 metadata。 |
| `scripts/build_glycosyl_donor_surrogates.py` | 沙盒/禁用 | 审核糖供体替代模板；当前模板均为 `placeholder_unmapped`，不得启用。 |
| `scripts/run_panel_ablation_experiments.py` | 活跃/待运行 | 对 PubChem 验证的黄酮苷 panel 运行 8/9 号虚拟桥接实验。 |
| `scripts/compare_panel_ablation_results.py` | 活跃/待运行 | 比较 panel 实验结果，并保持 bridge/ZINC/strict-trusted 分类边界。 |
| `scripts/build_strict_stock.py` | 活跃/支持 | 构建严格可购买库存并分离未验证的键。 |

## 活跃配置

| 路径 | 状态 | 作用 |
|---|---|---|
| `config/reaction_families.json` | 活跃 | 已审核黄酮反应族模板的来源。 |
| `config/sugar_closure_templates.json` | 活跃策略清单 | 记录糖闭环策略和未来禁用的模板族。 |
| `config/glycosyl_donor_surrogates.json` | 沙盒/禁用 | 未映射糖供体替代模板配置；全部 `active_expansion=false`，`validation_status=failed_map_retention`。 |
| `config/flavonoid_target_panel.csv` | 活跃 panel 配置 | PubChem 验证的 hesperidin、naringin、neohesperidin、narirutin、rutin、quercitrin 靶标结构。 |
| `config/ablation_*.yml` | 活跃实验 | ZINC/严格/可信/虚拟桥接测试的九个消融配置。 |
| `config/flavonoid.yml` | 旧版/兼容性 | 较旧的黄酮配置；除非明确刷新，否则不要用于当前结论。 |
| `config/hesperidin_optimized.yml` | 旧版/兼容性 | 较旧的优化橙皮苷配置；除非明确刷新，否则不要用于当前结论。 |
| `config/baseline.yml` | 旧版基线 | 较旧的基线配置。 |

## 活跃模板和库存产物

| 路径 | 状态 | 作用 |
|---|---|---|
| `templates/reaction_families/flavonoid_reaction_family_templates.hdf5` | 活跃生成 | MCTS 搜索库，仅包含活跃的反应族模板。 |
| `templates/reaction_families/flavonoid_reaction_family_templates_audit.csv` | 活跃审核 | 完整的活跃/非活跃模板审核表。 |
| `templates/stock_layers/strict_buyable_stock_inchikeys.txt` | 活跃库存 | 严格库存层。 |
| `templates/stock_layers/trusted_intermediate_stock_inchikeys.txt` | 活跃库存 | 可信中间层。 |
| `templates/stock_layers/virtual_bridge_stock_inchikeys.txt` | 活跃诊断库存 | 仅用于连通性诊断的虚拟桥接库存。 |
| `templates/stock_layers/sugar_gap_clusters.csv` | 活跃审核 | 路由间隙糖叶子的分类。 |
| `templates/stock_layers/sugar_bridge_stock.csv` | 活跃审核/生成库存 | 路由间隙衍生的糖桥接条目。 |
| `templates/stock_layers/sugar_bridge_evidence_review.csv` | 活跃审核 | 17 个 sugar bridge 条目的证据等级、保护态、连接候选和保守决策表。 |
| `templates/stock_layers/sugar_bridge_core_assignments.csv` | 活跃审核 | 保护态 sugar bridge 去乙酰化后的 sugar core 归一化分组。 |
| `templates/stock_layers/rutinosyl_chloride_structure_candidates.csv` | 活跃审核/沙盒 | 两个重建的 hexa-O-acetyl-rutinosyl chloride 异头候选；不作为库存输入。 |
| `templates/stock_layers/rutinosyl_chloride_bridge_mapping.csv` | 活跃审核/沙盒 | donor 候选折回 `C12H22O9` bridge skeleton 的连通性映射审计表。 |
| `templates/stock_layers/rutinose_donor_evidence_worklist.csv` | 活跃审核 | activated rutinose donor 证据工作列表；不作为库存输入。 |
| `templates/stock_layers/stock_layers_metadata.csv` | 活跃元数据 | 库存层证据和角色元数据。 |

## 当前可信输出

| 路径 | 状态 | 作用 |
|---|---|---|
| `outputs/ablation/ablation_report.md` | 活跃报告 | 人类可读的最新消融摘要。 |
| `outputs/ablation/ablation_report.json` | 活跃报告 | 结构化的最新消融摘要。 |
| `outputs/ablation/route_gap_worklist.csv` | 活跃工作列表 | 当前桥接层后的最新未解析叶子工作列表。 |
| `logs/reaction_family_template_audit.md` | 活跃审核 | 反应族模板审核。 |
| `logs/sugar_bridge_layer_audit.md` | 活跃审核 | 糖桥接层审核。 |
| `logs/sugar_bridge_evidence_review.md` | 活跃审核 | sugar bridge 证据等级和保护态分布。 |
| `logs/sugar_bridge_core_assignment.md` | 活跃审核 | sugar bridge 去乙酰化归一到 rutinose-like anomeric-deoxy bridge skeleton 的分组报告。 |
| `logs/rutinosyl_chloride_reconstruction_audit.md` | 活跃审核/沙盒 | hexa-O-acetyl-rutinosyl chloride 两个机器候选的 formula/乙酰/Cl/OH 验证报告。 |
| `logs/rutinosyl_chloride_bridge_mapping_audit.md` | 活跃审核/沙盒 | chloride donor 候选与 `C12H22O9` bridge artifact 的连通性映射报告。 |
| `logs/rutinose_donor_evidence_audit.md` | 活跃审核 | rutinose donor 文献/结构证据审计；记录一个命名 activated donor 文献命中和未升级决策。 |
| `logs/disaccharide_linkage_analysis.md` | 活跃审核 | rutinose/neohesperidose 连接与 InChIKey 区分。 |
| `logs/glycosyl_donor_surrogate_audit.md` | 沙盒审核 | donor surrogate 全部未映射、map retention 0、禁用状态的审计。 |
| `docs/execution_report.md` | 活跃记录 | 从虚拟糖桥到可审计供体层的执行总结。 |
| `docs/self_check_report.md` | 活跃记录 | AI 执行后的自检结果；需结合本 manifest 的人工复核口径解读。 |

## 旧版区域

| 路径 | 状态 | 备注 |
|---|---|---|
| `src/` | 旧版模板实验 | 早期提取/转换/诊断脚本。保留用于溯源，但未经审查不要用于当前结论。 |
| `templates/flavonoid_templates.*` | 旧版自定义库 | 较旧的模板集；当前 MCTS 消融使用反应族模板。 |
| `templates/flavonoid_structural_diversity/` | 旧版模板族 | 保留用于溯源的结构多样性模板。 |
| `templates/flavonoid_biosynthesis/` | 旧版模板族 | 保留用于溯源的生物合成模板。 |
| `templates/literature_curated/` | 证据/溯源 | 文献衍生的模板和证据表。`reaction_family_evidence_minimal.csv` 仍是相关证据。 |
| `archive/` | 归档输出/日志 | 较旧的路由输出和日期报告已移出活跃路径。 |

## 当前结果快照

糖桥接层和证据审计后：

| 实验 | 有效已解决 | strict/trusted | bridge-closed | ZINC baseline | 非虚拟 |
|---|---:|---:|---:|---:|---:|
| `baseline_zinc` | 10 | 0 | 0 | 10 | 10 |
| `flavonoid_zinc` | 3 | 0 | 0 | 3 | 3 |
| `custom_only_virtual_bridge` | 2 | 0 | 2 | 0 | 0 |
| `flavonoid_virtual_bridge` | 10 | 0 | 10 | 0 | 0 |
| `flavonoid_strict` | 0 | 0 | 0 | 0 | 0 |
| `flavonoid_trusted` | 0 | 0 | 0 | 0 | 0 |

解读：已审核的苷键切割模板达到了正确的糖/苷元断键。已解决的虚拟桥接路由验证了连通性，但不能证明严格可购买性。ZINC 路线单独作为数据库基线，不计入 strict/trusted solved。

## 当前证据审计状态

- `sugar_bridge_stock.csv` 中 17 个 route-gap 衍生条目没有自动提升到 strict/trusted，均仍只允许进入 `virtual_bridge`。
- 17 个条目中 16 个为乙酰化糖片段，1 个为未保护 bridge skeleton；该未保护条目 `UZIKLNYKVUKZQZ-IFLAJBTPSA-N` 的分子式是 `C12H22O9`，而 true rutinose 是 `C12H22O10`，因此它不是 free rutinose，而是缺少异头位氧/离去基信息的 rutinose-like anomeric-deoxy route-gap artifact。
- `UZIKLNYKVUKZQZ-IFLAJBTPSA-N` 与 metadata/PubChem 的 rutinose InChIKey 不同，不能再仅解释为还原端环式/链式平衡；当前正确解释是 bridge cleavage 产生了少一个氧的 sugar skeleton。它只能保留为连通性证据，不得作为 free sugar、donor seed 或 strict/trusted stock。
- 16 个乙酰化 sugar bridge 条目全部去乙酰化归一到 `UZIKLNYKVUKZQZ-IFLAJBTPSA-N`，因此视为同一 protected/anomeric-deoxy bridge artifact family，而不是 16 个独立证据候选，也不是 protected true-rutinose 证据；当前无需逐个手工审查。
- rutinose donor 证据工作列表已建立：当前最强直接命中是 `hexa-O-acetyl-beta-rutinosyl chloride` 文献候选（DOI: `10.1016/S0008-6215(00)84374-8`）。本地已从 true rutinose cyclic seed 重建两个 `C24H33ClO15` hexa-O-acetyl-rutinosyl chloride 异头候选，二者均为 6 个 O-acetyl、0 个残余 OH、1 个 Cl；但 beta assignment 尚未由 primary source 确认，因此不能进入 strict/trusted stock。
- 两个重建 chloride donor 候选去乙酰并将 Cl 归一为 H 后，均折回 `C12H22O9` bridge skeleton 的同一 InChIKey connectivity block；完整 stereochemical InChIKey 不同，因此当前只证明 donor/bridge 连通性对应，不证明 route-gap bridge 的精确立体结构身份。
- `rutinosyl trichloroacetimidate`、`rutinosyl bromide/acetobromorutinose` 的直接检索没有得到可用 rutinosyl donor 结构命中；这些 donor 类仍保留为 disabled placeholder，不能凭通用糖基化化学启用。
- `rutinose` 与 `neohesperidose` 已在 `stock_layers_metadata.csv` 中使用不同结构和不同 InChIKey：
  - rutinose: `OVVGHDNPYGTYIT-MBXIIVTHSA-N`
  - neohesperidose: `LTHOGZQBHZQCGR-CUHPIOEGSA-N`
- donor surrogate 模板全部为未映射占位：SMARTS 可解析，但 atom map 缺失、map retention 为 0、`active_expansion=false`。
- 黄酮苷 panel 结构来自 PubChem 并经 RDKit InChIKey 校验；尚未纳入当前 ablation 主结论。

## 下一步科学任务

1. 获取 primary paper 或等价结构来源，确认 `hexa-O-acetyl-beta-rutinosyl chloride` 对应两个重建 anomer 候选中的哪一个；未确认前不得声称 beta donor exact identity。
2. 在 beta/anomer 候选确认后，将 donor surrogate 从 `placeholder_unmapped` 重写为 atom-mapped 模板，并在 map-retention、结构身份和测试分子上通过审计后再考虑启用。
3. 将 `C12H22O9` bridge artifact 到 `C24H33ClO15` donor 候选的连通性映射继续保持为 sandbox 证据；除非完整立体结构被验证，否则不得作为 exact identity 证据。
4. 继续检索 rutinose donor 的供应商/文献全文证据；只有 exact structure、来源和用途都充分者才可进入 trusted/strict 候选。
5. 运行 PubChem 验证 panel 的 8/9 号虚拟桥接实验，评估 sugar bridge 对 naringin、neohesperidin、narirutin、rutin、quercitrin 的泛化。
6. 对 USPTO 保护糖 artifact 建立惩罚、归一化或过滤策略，避免保护态糖片段在报告中被误解为 donor 证据。
7. 将芳香黄酮苷叶子保持为人工审查目标，而非虚拟库存。

## 最亟待解决的问题（按优先级列）

1. **strict/trusted 仍然是 0 solved。**
   现在 `flavonoid_virtual_bridge = 10/10 solved`，但全部是 `bridge-closed`，不是 strict buyable solved。它证明路线连通性，但还不能证明真实可合成。

2. **17 个 sugar bridge 条目已完成保守分级和 bridge skeleton 归一化，但仍缺可升级证据。**
   现在确认 17 条全部归一到同一个 `C12H22O9` rutinose-like anomeric-deoxy bridge artifact，而不是 true rutinose `C12H22O10`。donor worklist 已找到一个命名 activated donor 文献线索，并生成两个 `C24H33ClO15` chloride donor 机器候选；但 beta/anomer assignment 和供应商证据仍未确认，不能仅凭 bridge-closed、糖芯识别、文献标题或本地重建升级。

3. **剩余 gap 已经转移到芳香黄酮糖苷片段。**
   sugar core 闭合后，未解决 leaf 从 53 降到 32，Top gap 变成 `manual_review_aromatic_glycoside`。这些不能再自动加 stock，否则会把“目标分子/保护态目标分子”误当起始物。

4. **真实 glycosyl donor 层仍未通过验证，但优先候选已收敛。**
   目前模板是：
   ```text
   flavonoid glycoside -> aglycone + neutral sugar fragment
   ```
   但真实化学需要：
   ```text
   aglycone phenol + activated sugar donor -> glycoside
   ```
   donor surrogate / trichloroacetimidate / bromide / thioglycoside 已有禁用占位配置，但全部未映射、map retention=0，不能进入 MCTS 结论。当前已完成 `hexa-O-acetyl-rutinosyl chloride` 两个 anomer 候选的机器重建；下一步是 primary-source 确认 beta 候选，然后再写 atom-mapped donor sandbox，而不是继续泛化增加 donor 类。

5. **rutinose / neohesperidose 已在 metadata 中区分，但 bridge 条目的具体连接仍多为 unknown。**
   两个标准糖记录已分开；route-gap 衍生的 17 个 bridge 条目仍多为 `rhamnosyl_position_unknown`，需要更细的连接识别。

6. **chalcone 模板仍然是 disabled。**
   它们目前不能通过 map-retention 和 test molecule 匹配。黄酮通用路线如果要走 chalcone/flavanone 环化，必须单独重写 chalcone sandbox。

7. **USPTO 保护基分支会制造大量糖保护态 artifact。**
   很多 `0.0 Unrecognized` 上游来自 USPTO 的保护/脱保护幻想路线。我们需要对保护糖片段做归一化或惩罚，否则搜索会偏向假保护态。

8. **bridge-closed solved 已单独分类，后续必须继续保持。**
   当前报告使用 `strict_trusted_solved`、`bridge_closed_connectivity`、`zinc_baseline_solved`、`unsolved` 四类，后续路线排序、结论文字、图表都要继续保持这个边界。

9. **跨黄酮测试集已配置但尚未运行成为主结论。**
   `config/flavonoid_target_panel.csv` 已使用 PubChem 结构并通过 InChIKey 校验；下一步需要运行 panel ablation 并人工审计结果。

10. **下一步真正的科学问题是“beta donor 精确结构确认 + atom-mapped donor sandbox”，不是继续扩 virtual stock。**
   virtual bridge 已经证明瓶颈在糖层。再盲目加 stock 会让结果变松；现在已有两个 `hexa-O-acetyl-rutinosyl chloride` 机器候选和 bridge connectivity 映射，下一步要确认哪一个是文献 beta donor，并把它变成可审计的映射模板测试结果，同时防止 `C12H22O9` bridge skeleton 被误读成 true rutinose。
