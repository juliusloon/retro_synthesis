# 第 3/4 件事执行计划：Panel 泛化与保护态 Sugar Artifact 控制

最后更新：2026-07-01

## 背景与边界

当前主线已经完成两件关键前置工作：

1. `hexa-O-acetyl-beta-rutinosyl chloride` 的 CW 候选已被人工确认为 beta donor stereochemical candidate。
2. CW beta donor 已进入 paired inactive atom-mapped sandbox；`active_expansion=false`，不进入 MCTS 主结论。

因此下一步不应继续扩展 virtual stock，也不应启用 donor sandbox。当前要解决的是：

- **第 3 件事**：运行 PubChem 验证的黄酮苷 panel 8/9 实验，评估现有 virtual bridge 对不同糖型和不同黄酮母核的泛化。
- **第 4 件事**：建立保护态 sugar artifact 的归一化、惩罚和报告口径，防止 USPTO 保护/脱保护幻想路线被误读为真实 donor 证据。

核心边界必须继续保持：

```text
strict/trusted solved
bridge-closed connectivity
ZINC baseline solved
unsolved / manual-review gaps
```

`bridge-closed` 只能证明连通性闭合，不能证明可购买性或真实 donor 可用性。

## 第 3 件事：PubChem Panel 8/9 泛化实验

### 目标

评估当前已审核 flavonoid templates + virtual bridge layer 在跨黄酮 panel 上是否稳定：

| 靶标 | 类型 | 预期糖型 | 目的 |
|---|---|---|---|
| hesperidin | flavanone glycoside | rutinose | 主线阳性对照 |
| narirutin | flavanone glycoside | rutinose | 同糖型、不同 aglycone |
| rutin | flavonol glycoside | rutinose | 同糖型、不同母核 |
| naringin | flavanone glycoside | neohesperidose | 糖连接差异测试 |
| neohesperidin | flavanone glycoside | neohesperidose | 糖连接差异测试 |
| quercitrin | flavonol glycoside | rhamnoside | 单糖负/边界测试 |

重点不是追求 strict/trusted solved，而是回答：

1. `flavonoid_virtual_bridge` 是否能在 rutinose / neohesperidose / rhamnoside 上合理闭合？
2. `custom_only_virtual_bridge` 是否也能闭合，还是依赖 USPTO 模板？
3. 未闭合 leaf 是否集中转移到 `manual_review_aromatic_glycoside` 或特定糖型？
4. 是否出现“目标分子/保护态目标分子被当成库存”的误判风险？

### 运行前检查

先确认 reaction family、stock layer、donor sandbox 审计处于最新状态：

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_reaction_family_templates.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_sugar_bridge_layer.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/audit_sugar_bridge_evidence.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/group_sugar_bridge_cores.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_glycosyl_donor_surrogates.py
```

检查点：

- `config/flavonoid_target_panel.csv` 仍为 `verified_pubchem`。
- `scripts/run_panel_ablation_experiments.py` 会将 panel SMILES canonicalize 后再喂给 AiZynthFinder。
- `config/glycosyl_donor_surrogates.json` 中 paired donor sandbox 仍为 `active_expansion=false`。

### 正式运行

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_panel_ablation_experiments.py --experiments 8,9
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_panel_ablation_results.py
```

预期输出：

```text
outputs/panel_ablation/panel_ablation_summary.json
outputs/panel_ablation/panel_ablation_report.md
outputs/panel_ablation/panel_ablation_report.json
outputs/panel_ablation/<target>/flavonoid_virtual_bridge.json
outputs/panel_ablation/<target>/custom_only_virtual_bridge.json
```

### 分析项目

对每个 target / experiment 记录：

| 指标 | 解释 | 必须保留的判断 |
|---|---|---|
| `n_solved` | AiZynth 原生 solved 数 | 不能单独作为科学结论 |
| `n_bridge_closed` | 使用 virtual bridge 闭合的 solved | 只算 connectivity |
| `n_non_virtual` | 不依赖 virtual bridge 的 solved | 才可能进入更严格讨论 |
| sugar family | rutinose / neohesperidose / rhamnoside | 判断泛化范围 |
| leaf InChIKey 分布 | 未解决叶子 | 判断 gap 是否转移 |
| route validity class | strict / bridge / ZINC / unsolved | 保持主报告分类边界 |

需要额外人工检查的情况：

- `rhamnoside` 或 `neohesperidose` 被错误闭合到 rutinose-like bridge。
- route leaf 是完整靶标、保护态靶标或高度相似的芳香黄酮苷片段。
- USPTO 分支产生多乙酰化 sugar fragment，并被报告成 donor 证据。
- `custom_only_virtual_bridge` 与 `flavonoid_virtual_bridge` 结论差异很大，说明 USPTO 主导了闭合。

### 验收标准

第 3 件事完成的标准：

- `panel_ablation_report.md/json` 已生成。
- 每个 target 的 8/9 实验都有结果或明确失败原因。
- 报告中显式区分：
  - `bridge_closed_connectivity`
  - `strict_trusted_solved`
  - `zinc_baseline_solved`
  - `unsolved`
- 对 rutinose / neohesperidose / rhamnoside 三类糖型分别给出一句结论。
- 若发现异常闭合，生成一个后续 worklist，而不是直接升级 stock。

### 不做的事

- 不启用 donor sandbox。
- 不把 `virtual_bridge` 结果写成 strict/trusted。
- 不因为 panel solved 数增加就宣传真实可合成。
- 不把 `C12H22O9` bridge skeleton 当作 true rutinose。

## 第 4 件事：保护态 Sugar Artifact 控制

### 目标

把 USPTO 保护/脱保护分支产生的保护态 sugar artifact 从“候选 donor 证据”中剥离出来，只作为 route-gap / search-bias 诊断。

当前已知风险：

- `0.0 Unrecognized` 上游经常产生乙酰化 sugar fragment。
- 17 个 sugar bridge 条目中 16 个乙酰化条目已归一到同一个 `C12H22O9` bridge artifact family。
- 这些条目不是 16 个独立 donor 候选，也不是 protected true-rutinose 证据。

### 工作包 A：Artifact 分类规则

在 route gap / sugar bridge 审计中统一以下标签：

| 标签 | 判据 | 报告解释 |
|---|---|---|
| `anomeric_deoxy_bridge_artifact` | 分子式为 `C12H22O9` 或归一后为该 skeleton | 连通性 artifact，不是 free sugar |
| `protected_bridge_artifact` | acetyl/silyl 等保护基数 > 0，且去保护后回到 bridge artifact | 保护态 route-gap artifact |
| `aromatic_glycoside_manual_review` | 仍含黄酮/芳香苷元片段 | 人工审查，不可自动加 stock |
| `candidate_real_donor` | 有离去基、exact donor 结构、来源证据 | 进入 donor evidence worklist，不直接入 stock |

### 工作包 B：归一化策略

对 sugar-like leaf 做两层归一化：

1. **保护基归一化**：去 acetyl / silyl 等常见保护基，记录原保护基数。
2. **bridge skeleton 归一化**：检查是否回到 `UZIKLNYKVUKZQZ-IFLAJBTPSA-N` connectivity block。

输出建议：

```text
templates/stock_layers/protected_sugar_artifact_review.csv
logs/protected_sugar_artifact_audit.md
```

建议字段：

```text
name
smiles
inchikey
source_experiment
upstream_family
acetyl_count
silyl_count
normalized_smiles
normalized_inchikey
normalized_family
artifact_class
allowed_report_role
stock_decision
notes
```

### 工作包 C：报告惩罚/过滤口径

比较报告中新增或强化三个统计：

| 统计 | 目的 |
|---|---|
| `contains_protected_sugar_artifact` | 标记路线是否经过保护态糖 artifact |
| `protected_artifact_leaf_count` | 量化 USPTO 保护态偏置 |
| `manual_review_aromatic_glycoside_count` | 防止目标样片被误当起始物 |

路线排序建议：

1. strict/trusted non-artifact
2. bridge-closed without protected artifact
3. bridge-closed with protected artifact
4. unsolved / manual review

保护态 artifact 不应让路线进入更高证据等级；它只能作为警告或惩罚项。

### 工作包 D：负面规则

以下情况必须明确拒绝升级：

- 只有 protected sugar fragment，没有 exact leaving group。
- 去保护后回到 `C12H22O9` bridge skeleton。
- 只来自 USPTO 保护/脱保护幻想路线，没有文献或供应商结构证据。
- 芳香黄酮苷片段仍像目标分子或保护态目标分子。
- 只凭 common glycosylation chemistry 推断 bromide / trichloroacetimidate / thioglycoside donor。

### 验收标准

第 4 件事完成的标准：

- 保护态 sugar artifact 有独立 audit 表和日志。
- 16 个乙酰化 sugar bridge 不再在报告中表现为 16 个 donor 候选。
- 报告中可以区分：
  - bridge skeleton artifact
  - protected bridge artifact
  - aromatic glycoside manual review
  - real donor evidence candidate
- route ranking 或 summary 中显式惩罚/标注 protected artifact。
- 没有任何 protected artifact 被提升到 strict/trusted stock。

## 建议执行顺序

1. 跑 panel 8/9，并生成 `panel_ablation_report.md/json`。
2. 人工读 panel report，记录 sugar family 泛化结论。
3. 从 panel route leaves 和现有 `route_gap_worklist.csv` 合并提取 artifact 候选。
4. 新建 protected sugar artifact 审计脚本和输出表。
5. 修改 compare 脚本，让 protected artifact 成为报告中的显式惩罚/警告。
6. 更新 `current_pipeline_manifest.md`，把 panel 结果和 artifact 控制结论写入活跃状态。

## 风险与决策点

| 风险 | 处理 |
|---|---|
| panel solved 全靠 virtual bridge | 保留为 connectivity，不写 strict |
| neohesperidose 被 rutinose bridge 错闭合 | 标为糖型错误，需要 bridge layer 拆分 |
| protected sugar 数量很多 | 不扩 stock，先归一化和惩罚 |
| donor sandbox 看起来能闭合路线 | 仍保持 inactive，除非 production gating 完成 |
| aromatic glycoside leaf 很像目标 | manual review，不自动库存化 |

## 最终交付物

第 3 件事交付：

```text
outputs/panel_ablation/panel_ablation_report.md
outputs/panel_ablation/panel_ablation_report.json
docs/current_pipeline_manifest.md 的 panel 结果更新
```

第 4 件事交付：

```text
templates/stock_layers/protected_sugar_artifact_review.csv
logs/protected_sugar_artifact_audit.md
compare/report 中的 protected artifact 惩罚或警告统计
docs/current_pipeline_manifest.md 的 artifact 控制结论更新
```

最终一句话标准：

```text
Panel 结果用于评估 bridge 泛化；protected sugar artifact 用于诊断搜索偏置。
二者都不能绕过 strict/trusted 证据门槛。
```
