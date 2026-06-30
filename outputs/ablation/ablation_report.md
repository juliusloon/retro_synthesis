# 消融实验报告

生成日期: 2026-06-30 19:59

## 评价标准

- **AiZynth solved**: AiZynthFinder 原生 `is_solved` 字段
- **map-valid solved**: AiZynth solved 且所有反应步骤 atom-map retention ≥ 0.8
- **有效 solved**: map-valid solved（即真正通过守恒审计的路线）
- **bridge-closed solved**: 有效 solved 且至少一个叶子节点来自 `virtual_bridge`，只用于连通性验证，不等同于真实 strict solved
- **non-virtual effective solved**: 有效 solved 且不依赖 `virtual_bridge`
- **stock 层级 solved**: 路线中所有叶子节点均属于该 stock 层级
- **关键中间体到达**: 路线中包含到达 aglycone/chalcone/sugar donor 的步骤

## 汇总

| 实验 | 总路线 | AiZynth | map-valid | 有效 solved | bridge-closed | non-virtual | 最短步 |
|---|---:|---:|---:|---:|---:|---:|---:|
| baseline_strict | 12 | 0 | 0 | 0 | 0 | 0 | - |
| baseline_zinc | 10 | 10 | 10 | 10 | 0 | 10 | 1 |
| custom_only_strict | 5 | 0 | 0 | 0 | 0 | 0 | - |
| custom_only_trusted | 4 | 0 | 0 | 0 | 0 | 0 | - |
| custom_only_virtual_bridge | 4 | 2 | 2 | 2 | 2 | 0 | 1 |
| flavonoid_strict | 5 | 0 | 0 | 0 | 0 | 0 | - |
| flavonoid_trusted | 7 | 0 | 0 | 0 | 0 | 0 | - |
| flavonoid_virtual_bridge | 10 | 10 | 10 | 10 | 10 | 0 | 1 |
| flavonoid_zinc | 7 | 3 | 3 | 3 | 0 | 3 | 1 |

## Stock 层级分析

| 实验 | strict_buyable | trusted_intermediate | virtual_bridge | unknown |
|---|---:|---:|---:|---:|
| baseline_strict | 0 | 0 | 0 | 0 |
| baseline_zinc | 2 | 0 | 0 | 10 |
| custom_only_strict | 0 | 0 | 0 | 0 |
| custom_only_trusted | 0 | 0 | 0 | 0 |
| custom_only_virtual_bridge | 1 | 2 | 2 | 0 |
| flavonoid_strict | 0 | 0 | 0 | 0 |
| flavonoid_trusted | 0 | 0 | 0 | 0 |
| flavonoid_virtual_bridge | 2 | 9 | 10 | 0 |
| flavonoid_zinc | 0 | 0 | 0 | 3 |

## 自定义模板家族统计 (all routes)

| 实验 | 自定义模板家族 | 使用次数 |
|---|---|---:|
| baseline_strict | - | 0 |
| baseline_zinc | - | 0 |
| custom_only_strict | aryl_methyl_ether_cleavage | 3 |
| custom_only_strict | o_glycoside_cleavage_pyranose | 3 |
| custom_only_trusted | aryl_methyl_ether_cleavage | 2 |
| custom_only_trusted | o_glycoside_cleavage_pyranose | 2 |
| custom_only_virtual_bridge | o_glycoside_cleavage_pyranose | 2 |
| custom_only_virtual_bridge | aryl_methyl_ether_cleavage | 2 |
| flavonoid_strict | o_glycoside_cleavage_pyranose | 5 |
| flavonoid_strict | aryl_methyl_ether_cleavage | 5 |
| flavonoid_trusted | aryl_methyl_ether_cleavage | 12 |
| flavonoid_trusted | o_glycoside_cleavage_pyranose | 7 |
| flavonoid_trusted | aryl_O_methylation | 2 |
| flavonoid_virtual_bridge | o_glycoside_cleavage_pyranose | 10 |
| flavonoid_zinc | o_glycoside_cleavage_pyranose | 4 |

## 关键中间体到达统计

| 实验 | aglycone (solved/all) | chalcone (solved/all) | sugar donor (solved/all) |
|---|---:|---:|---:|
| baseline_strict | 0/0 | 0/0 | 0/0 |
| baseline_zinc | 0/0 | 0/0 | 0/0 |
| custom_only_strict | 0/3 | 0/0 | 0/0 |
| custom_only_trusted | 0/2 | 0/0 | 0/0 |
| custom_only_virtual_bridge | 2/2 | 0/0 | 0/0 |
| flavonoid_strict | 0/5 | 0/0 | 0/0 |
| flavonoid_trusted | 0/7 | 0/0 | 0/0 |
| flavonoid_virtual_bridge | 10/10 | 0/0 | 0/0 |
| flavonoid_zinc | 0/4 | 0/0 | 0/0 |

## 路由有效性分类统计

| 实验 | strict_trusted_solved | bridge_closed_connectivity | zinc_baseline_solved | unsolved |
|---|---:|---:|---:|---:|
| baseline_strict | 0 | 0 | 0 | 12 |
| baseline_zinc | 0 | 0 | 10 | 0 |
| custom_only_strict | 0 | 0 | 0 | 5 |
| custom_only_trusted | 0 | 0 | 0 | 4 |
| custom_only_virtual_bridge | 0 | 2 | 0 | 2 |
| flavonoid_strict | 0 | 0 | 0 | 5 |
| flavonoid_trusted | 0 | 0 | 0 | 7 |
| flavonoid_virtual_bridge | 0 | 10 | 0 | 0 |
| flavonoid_zinc | 0 | 0 | 3 | 4 |

## 虚拟桥接使用统计

| 实验 | uses_virtual_bridge | uses_sugar_gap_bridge | contains_protected_sugar_artifact |
|---|---:|---:|---:|
| baseline_strict | 0 | 0 | 0 |
| baseline_zinc | 0 | 0 | 0 |
| custom_only_strict | 3 | 3 | 0 |
| custom_only_trusted | 2 | 2 | 0 |
| custom_only_virtual_bridge | 2 | 2 | 0 |
| flavonoid_strict | 0 | 0 | 0 |
| flavonoid_trusted | 6 | 6 | 6 |
| flavonoid_virtual_bridge | 10 | 10 | 7 |
| flavonoid_zinc | 4 | 4 | 4 |

## 结论

- Stock 影响 (Baseline): ZINC 有效=10, Strict 有效=0
- Custom template 影响 (Strict): Baseline 有效=0, +Custom 有效=0
- 仅 Custom template: AiZynth=0, 有效=0
- Sugar bridge (Custom only + virtual): 有效=2, bridge-closed=2, non-virtual=0
- Sugar bridge (USPTO+Custom + virtual): 有效=10, bridge-closed=10, non-virtual=0
