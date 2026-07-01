# 保护态 Sugar Artifact 审计报告

生成日期: 2026-07-01 16:57

## 概述

- 总 artifact 候选: 18
- 来自 sugar bridge evidence review: 17
- 来自 route gap worklist (sugar bridge): 1

## 分类统计

| Artifact 类 | 数量 | 报告角色 | 库存决策 |
|---|---:|---|---|
| anomeric_deoxy_bridge_artifact | 1 | connectivity_evidence_only | keep_virtual_bridge_only |
| protected_bridge_artifact | 17 | search_bias_diagnostic | keep_virtual_bridge_only |

## 保护基统计

- 含 O-acetyl 保护基: 17
- 含 silyl 保护基: 0
- 归一化到 bridge skeleton (C12H22O9): 18

## 负面规则验证

以下规则已自动检查：

1. **只有 protected sugar fragment，没有 exact leaving group** → 分为 `protected_bridge_artifact`，不升级
2. **去保护后回到 C12H22O9 bridge skeleton** → 分为 `protected_bridge_artifact` 或 `anomeric_deoxy_bridge_artifact`
3. **只来自 USPTO 保护/脱保护幻想路线** → `upstream_family` 标记来源
4. **芳香黄酮苷片段仍像目标分子** → 分为 `aromatic_glycoside_manual_review`
5. **只凭 common glycosylation chemistry 推断 donor** → 无 candidate_real_donor 条目

## anomeric_deoxy_bridge_artifact

| 名称 | InChIKey | Formula | Acetyl | 来源 | 库存决策 |
|---|---|---|---:|---|---|
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `UZIKLNYKVUKZQZ-IFLAJ...` | C12H22O9 | 0 | sugar_bridge_evidence_review | keep_virtual_bridge_only |

## protected_bridge_artifact

| 名称 | InChIKey | Formula | Acetyl | 来源 | 库存决策 |
|---|---|---|---:|---|---|
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `SPCUBGINFCWHNO-ZGYAB...` | C14H24O10 | 1 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `ASGYRTWPIYFTQG-UHBWO...` | C22H32O14 | 5 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `BSESJCOUYFELHO-DFRUX...` | C16H26O11 | 2 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `COAOCSAWRKFHQU-OYESW...` | C16H26O11 | 2 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `HMSBGGCYEHLLBH-NXBYO...` | C22H32O14 | 5 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `NDAQNSFGSPFDIX-JEUWS...` | C22H32O14 | 5 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `NTQPAWBRYIOBFQ-ZGYAB...` | C14H24O10 | 1 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `RAFSHHVKLARFJP-YOKAV...` | C22H32O14 | 5 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `SJFVSFLVJAZSFI-NXBYO...` | C22H32O14 | 5 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `TXZOTBUOBZIVMS-ZGYAB...` | C14H24O10 | 1 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `VQFCRWDUJQUJMR-MNENM...` | C22H32O14 | 5 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `MZLKGISSOYWIMT-ZFOBK...` | C14H24O10 | 1 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `OFZYAYMKZLYPAR-XZTLZ...` | C14H24O10 | 1 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `UDYLCPHGMFALPK-ZGYAB...` | C14H24O10 | 1 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `UNBOZMFZFLXXOQ-LIEZG...` | C16H26O11 | 2 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ | `VOAJQPZXTFIBIC-GRGQJ...` | C16H26O11 | 2 | sugar_bridge_evidence_review | keep_virtual_bridge_only |
| worklist_sugar_bridge_IHXLRXWU | `IHXLRXWUIMEDJT-DEFHS...` | C18H28O12 | 3 | uspto_custom_zinc | keep_virtual_bridge_only |

## 路线排序建议

1. strict/trusted non-artifact
2. bridge-closed without protected artifact
3. bridge-closed with protected artifact
4. unsolved / manual review

保护态 artifact 不应让路线进入更高证据等级；它只能作为警告或惩罚项。

## 结论

- 16 个乙酰化 sugar bridge 条目归一到 UZIK bridge artifact family（C12H22O9）
- 2 个 worklist 条目归一化后不是 UZIK bridge skeleton，而是 `CBHXWLSZNTXSTO-IBLCKKAJSA-N`（C18H34O9），分为 `protected_nonbridge_sugar_artifact`
- 不是独立 donor 候选，也不是 protected true-rutinose 证据
- 所有 protected artifact 保持在 virtual_bridge，不进入 strict/trusted
- aromatic glycoside 条目保持为人工审查目标
