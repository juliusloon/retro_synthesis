# 消融实验报告

生成日期: 2026-07-01 16:14

## 评价框架

当前 strict/trusted 库规模不足，不适合作为路线发现阶段的主成功标准；本阶段以 ZINC baseline 评估通用库存闭合能力，以 virtual bridge 诊断糖层连通性，以 strict/trusted 标记高置信库存子集。

| 层级 | 角色 |
|---|---|
| ZINC baseline | 主搜索基准，回答"大库存下路线能不能闭合" |
| virtual_bridge | 诊断糖层连通性瓶颈 |
| strict/trusted | 保守证据层，用来标注哪些叶子证据更硬 |
| donor sandbox | 未来生产级糖供体模板验证，不进主结论 |

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
| custom_only_full_stock | 4 | 3 | 3 | 3 | 1 | 2 | 1 |
| custom_only_zinc | 4 | 0 | 0 | 0 | 0 | 0 | - |
| uspto_custom_zinc | 12 | 3 | 3 | 3 | 0 | 3 | 1 |
| uspto_custom_zinc_strict | 6 | 4 | 4 | 4 | 0 | 4 | 1 |
| uspto_custom_zinc_trusted | 8 | 4 | 4 | 4 | 0 | 4 | 1 |
| uspto_custom_zinc_vbridge | 13 | 13 | 13 | 13 | 9 | 4 | 1 |
| uspto_rb_zinc | 5 | 5 | 5 | 5 | 0 | 5 | 1 |
| uspto_zinc | 10 | 10 | 10 | 10 | 0 | 10 | 1 |

## Stock 层级分析

| 实验 | strict_buyable | trusted_intermediate | virtual_bridge | unknown |
|---|---:|---:|---:|---:|
| custom_only_full_stock | 2 | 1 | 1 | 2 |
| custom_only_zinc | 0 | 0 | 0 | 0 |
| uspto_custom_zinc | 0 | 0 | 0 | 3 |
| uspto_custom_zinc_strict | 1 | 0 | 0 | 4 |
| uspto_custom_zinc_trusted | 1 | 0 | 0 | 4 |
| uspto_custom_zinc_vbridge | 2 | 8 | 9 | 4 |
| uspto_rb_zinc | 2 | 0 | 0 | 5 |
| uspto_zinc | 2 | 0 | 0 | 10 |

## 自定义模板家族统计 (all routes)

| 实验 | 自定义模板家族 | 使用次数 |
|---|---|---:|
| custom_only_full_stock | o_glycoside_cleavage_pyranose | 1 |
| custom_only_full_stock | aryl_methyl_ether_cleavage | 1 |
| custom_only_full_stock | aryl_O_methylation | 1 |
| custom_only_zinc | o_glycoside_cleavage_pyranose | 1 |
| custom_only_zinc | aryl_methyl_ether_cleavage | 1 |
| custom_only_zinc | aryl_O_methylation | 1 |
| uspto_custom_zinc | o_glycoside_cleavage_pyranose | 9 |
| uspto_custom_zinc_strict | o_glycoside_cleavage_pyranose | 2 |
| uspto_custom_zinc_strict | aryl_methyl_ether_cleavage | 2 |
| uspto_custom_zinc_trusted | o_glycoside_cleavage_pyranose | 4 |
| uspto_custom_zinc_trusted | aryl_methyl_ether_cleavage | 2 |
| uspto_custom_zinc_trusted | aryl_O_methylation | 2 |
| uspto_custom_zinc_vbridge | o_glycoside_cleavage_pyranose | 9 |
| uspto_rb_zinc | - | 0 |
| uspto_zinc | - | 0 |

## 关键中间体到达统计

| 实验 | aglycone (solved/all) | chalcone (solved/all) | sugar donor (solved/all) |
|---|---:|---:|---:|
| custom_only_full_stock | 1/1 | 0/0 | 0/0 |
| custom_only_zinc | 0/1 | 0/0 | 0/0 |
| uspto_custom_zinc | 0/9 | 0/0 | 0/0 |
| uspto_custom_zinc_strict | 0/2 | 0/0 | 0/0 |
| uspto_custom_zinc_trusted | 0/4 | 0/0 | 0/0 |
| uspto_custom_zinc_vbridge | 9/9 | 0/0 | 0/0 |
| uspto_rb_zinc | 0/0 | 0/0 | 0/0 |
| uspto_zinc | 0/0 | 0/0 | 0/0 |

## 路由有效性分类统计

| 实验 | strict_trusted_solved | bridge_closed_connectivity | zinc_baseline_solved | unsolved |
|---|---:|---:|---:|---:|
| custom_only_full_stock | 0 | 1 | 2 | 1 |
| custom_only_zinc | 0 | 0 | 0 | 4 |
| uspto_custom_zinc | 0 | 0 | 3 | 9 |
| uspto_custom_zinc_strict | 0 | 0 | 4 | 2 |
| uspto_custom_zinc_trusted | 0 | 0 | 4 | 4 |
| uspto_custom_zinc_vbridge | 0 | 9 | 4 | 0 |
| uspto_rb_zinc | 0 | 0 | 5 | 0 |
| uspto_zinc | 0 | 0 | 10 | 0 |

## 虚拟桥接使用统计

| 实验 | uses_virtual_bridge | uses_sugar_gap_bridge | contains_protected_sugar_artifact |
|---|---:|---:|---:|
| custom_only_full_stock | 1 | 1 | 0 |
| custom_only_zinc | 1 | 1 | 0 |
| uspto_custom_zinc | 8 | 8 | 8 |
| uspto_custom_zinc_strict | 2 | 2 | 2 |
| uspto_custom_zinc_trusted | 4 | 4 | 4 |
| uspto_custom_zinc_vbridge | 9 | 9 | 7 |
| uspto_rb_zinc | 0 | 0 | 0 |
| uspto_zinc | 0 | 0 | 0 |

## 保护态 Sugar Artifact 惩罚统计

保护态 sugar artifact 不应让路线进入更高证据等级；它只能作为警告或惩罚项。

| 实验 | 含保护态 artifact 路线数 | 含保护态 artifact solved 路线数 | 含芳香苷叶子路线数 | 保护态叶子总数 |
|---|---:|---:|---:|---:|
| custom_only_full_stock | 0 | 0 | 0 | 1 |
| custom_only_zinc | 0 | 0 | 0 | 1 |
| uspto_custom_zinc | 8 | 0 | 0 | 8 |
| uspto_custom_zinc_strict | 2 | 0 | 0 | 2 |
| uspto_custom_zinc_trusted | 4 | 0 | 0 | 4 |
| uspto_custom_zinc_vbridge | 7 | 7 | 0 | 9 |
| uspto_rb_zinc | 0 | 0 | 0 | 0 |
| uspto_zinc | 0 | 0 | 0 | 0 |

## 结论

- 结果集模式: `A1-C2`。
- 当前报告按 A1-C2 新实验框架解释；legacy 输出若同时存在，只保留在汇总表中，不用于主增益结论。
- 模板主比较候选: A1 effective=10, A3 effective=3。
- B3 full-stock/virtual-bridge 诊断: effective=13, bridge-closed=9, non-virtual=4。
- Custom-only 诊断: C1 effective=0, C2 effective=3。

## 增益分析

| 实验 | vs | 标签 | solved | 基线 solved | 增益 |
|---|---|---|---:|---:|---:|
| A2 | A1 | RingBreaker 增益 | 5 | 10 | -5 |
| A3 | A2 | Custom 模板增益 | 3 | 5 | -2 |
| B1 | A3 | strict 增益 | 4 | 3 | +1 |
| B2 | B1 | trusted 增益 | 4 | 4 | +0 |
| B3 | B2 | virtual bridge 增益 | 13 | 4 | +9 |
| C1 | A1 | 模板独立性（Custom only vs USPTO） | 0 | 10 | -10 |
| C2 | C1 | 全库存增益（Custom only） | 3 | 0 | +3 |
