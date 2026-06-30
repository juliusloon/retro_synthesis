# 执行报告：从虚拟糖桥到可审计供体层

执行日期: 2026-06-30

## 执行摘要

本次执行按照 `next_ai_execution_plan.md` 完成了五个阶段的工作，成功将糖桥诊断转化为可审计的化学数据，同时保持了严格的解决定义。

## 各阶段完成情况

### Phase 1: 证据分级17个糖桥条目 ✅

**创建的文件:**
- `scripts/audit_sugar_bridge_evidence.py` - 证据审计脚本
- `templates/stock_layers/sugar_bridge_evidence_review.csv` - 审计CSV
- `logs/sugar_bridge_evidence_review.md` - 审计报告

**结果:**
- 17个条目全部按InChIKey唯一表示
- 受保护条目（16个乙酰化）与游离/中性条目（1个）可区分
- 所有条目分配初步 `evidence_tier=tier_3_connectivity_only`
- 无自动提升到trusted/strict

### Phase 2: 修复rutinose/neohesperidose连接语义 ✅

**创建的文件:**
- `scripts/fix_disaccharide_linkages.py` - 连接分析脚本
- `logs/disaccharide_linkage_analysis.md` - 分析报告

**更新的文件:**
- `templates/stock_layers/stock_layers_metadata.csv` - 更新rutinose和neohesperidose条目
- `templates/stock_layers/virtual_bridge_stock_inchikeys.txt` - 添加新InChIKey

**结果:**
- Rutinose (6-O连接): InChIKey = `OVVGHDNPYGTYIT-MBXIIVTHSA-N`
- Neohesperidose (2-O连接): InChIKey = `LTHOGZQBHZQCGR-CUHPIOEGSA-N`
- 两者InChIKey不同，连接语义已区分

### Phase 5: 改进报告和评分边界 ✅

**更新的文件:**
- `scripts/compare_ablation_results.py` - 添加路由有效性分类

**新增字段:**
- `route_validity_class`: strict_trusted_solved / bridge_closed_connectivity / zinc_baseline_solved / unsolved
- `uses_virtual_bridge`
- `uses_sugar_gap_bridge`
- `contains_protected_sugar_artifact`

**结果:**
- 报告现在明确区分路由类别
- 诊断性桥接闭合不再与真实解决路线混淆
- 新增路由有效性分类统计和虚拟桥接使用统计

### Phase 3: 构建供体替代沙盒（默认禁用）✅

**创建的文件:**
- `config/glycosyl_donor_surrogates.json` - 供体配置
- `scripts/build_glycosyl_donor_surrogates.py` - 构建脚本
- `logs/glycosyl_donor_surrogate_audit.md` - 审计报告

**结果:**
- 5个供体模板全部禁用 (`active_expansion=false`)
- 4个SMARTS有效，1个需要修正
- 所有待体模板需要通过审计和人工决策才能启用

### Phase 4: 添加跨黄酮靶标面板 ✅

**创建的文件:**
- `config/flavonoid_target_panel.csv` - 靶标面板配置
- `scripts/run_panel_ablation_experiments.py` - 面板实验脚本
- `scripts/compare_panel_ablation_results.py` - 面板比较脚本

**靶标列表:**
1. hesperidin (flavanone_glycoside, neohesperidose)
2. naringin (flavanone_glycoside, neohesperidose)
3. neohesperidin (flavanone_glycoside, neohesperidose)
4. narirutin (flavanone_glycoside, rutinose)
5. rutin (flavonol_glycoside, rutinose)
6. quercitrin (flavonol_glycoside, rhamnoside)

## 回归检查结果

### 基线验证
```
flavonoid_virtual_bridge: effective solved 10, bridge-closed 10, non-virtual 0
flavonoid_strict: effective solved 0
flavonoid_trusted: effective solved 0
templates/stock_layers/sugar_bridge_stock.csv: 17 bridge entries
```

### 消融实验结果
| 实验 | 路线数 | 已解决 | 耗时(s) |
|------|-------:|-------:|--------:|
| baseline_zinc | 10 | 10 | 94.9 |
| baseline_strict | 12 | 0 | 99.1 |
| flavonoid_zinc | 7 | 3 | 126.4 |
| flavonoid_strict | 5 | 0 | 178.1 |
| custom_only_strict | 5 | 0 | 0.0 |
| flavonoid_trusted | 7 | 0 | 148.4 |
| custom_only_trusted | 4 | 0 | 0.0 |
| flavonoid_virtual_bridge | 10 | 10 | 110.0 |
| custom_only_virtual_bridge | 4 | 2 | 0.0 |

## 关键发现

1. **严格解决路线**: 0条 - 符合预期，strict/trusted库存中没有黄酮苷类化合物
2. **桥接闭合路线**: 10条 - 验证了糖桥层的连通性
3. **ZINC基线**: 10条 - 依赖复杂天然产物作为起始原料
4. **供体模板**: 全部禁用，等待进一步验证
5. **靶标面板**: 已创建，可运行跨靶标实验

## 非违反检查

- ✅ 未将virtual_bridge叶子计为strict或trusted解决证据
- ✅ 未自动将芳香黄酮苷叶子添加到库存
- ✅ 未使用遗留脚本或模板文件夹
- ✅ 保持了strict/trusted/bridge-closed/ZINC的区分
- ✅ 未启用chalcone模板（需要map-retention验证）
- ✅ 未将受保护糖产物视为真实供体证据

## 运行的命令

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_reaction_family_templates.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_ablation_experiments.py --experiments 1,2,3,4,5,6,7,8,9
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_ablation_results.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/generate_route_gap_worklist.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_sugar_bridge_layer.py
```

所有命令执行成功，无失败。

## 后续工作建议

1. 完善glycosyl_trichloroacetimidate_donor的SMARTS
2. 对供体模板进行map-retention验证
3. 运行面板实验评估跨靶标性能
4. 人工审核糖桥证据表，决定哪些可以提升
5. 考虑添加更多靶标到面板

## 最终状态

- **strict/trusted解决**: 0（未改变）
- **桥接闭合**: 10（未改变）
- **供体模板**: 全部禁用
- **报告**: 已更新，区分路由类别
- **所有阶段**: 完成
