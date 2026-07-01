# 跨黄酮靶标面板消融实验报告

生成日期: 2026-07-01 14:55

## 评价框架

当前 strict/trusted 库规模不足，不适合作为路线发现阶段的主成功标准；本阶段以 ZINC baseline 评估通用库存闭合能力，以 virtual bridge 诊断糖层连通性，以 strict/trusted 标记高置信库存子集。

| 层级 | 角色 |
|---|---|
| ZINC baseline | 主搜索基准，回答"大库存下路线能不能闭合" |
| virtual_bridge | 诊断糖层连通性瓶颈 |
| strict/trusted | 保守证据层，用来标注哪些叶子证据更硬 |
| donor sandbox | 未来生产级糖供体模板验证，不进主结论 |

## 靶标列表

| 靶标 | 类别 | 预期糖家族 |
|------|------|------------|
| hesperidin | flavanone_glycoside | rutinose |
| naringin | flavanone_glycoside | neohesperidose |
| narirutin | flavanone_glycoside | rutinose |
| neohesperidin | flavanone_glycoside | neohesperidose |
| quercitrin | flavonol_glycoside | rhamnoside |
| rutin | flavonol_glycoside | rutinose |

## 详细结果

### hesperidin

- 类别: flavanone_glycoside
- 预期糖家族: rutinose

| 实验 | 总路线 | ZINC solved | bridge-closed | strict/trusted |
|------|-------:|------------:|--------------:|---------------:|
| custom_only_virtual_bridge | 5 | 0 | 0 | 0 |
| flavonoid_virtual_bridge | 11 | 0 | 11 | 0 |
| flavonoid_zinc | 6 | 4 | 0 | 0 |

### naringin

- 类别: flavanone_glycoside
- 预期糖家族: neohesperidose

| 实验 | 总路线 | ZINC solved | bridge-closed | strict/trusted |
|------|-------:|------------:|--------------:|---------------:|
| custom_only_virtual_bridge | 2 | 0 | 0 | 0 |
| flavonoid_virtual_bridge | 12 | 0 | 0 | 0 |
| flavonoid_zinc | 13 | 13 | 0 | 0 |

### narirutin

- 类别: flavanone_glycoside
- 预期糖家族: rutinose

| 实验 | 总路线 | ZINC solved | bridge-closed | strict/trusted |
|------|-------:|------------:|--------------:|---------------:|
| custom_only_virtual_bridge | 2 | 0 | 0 | 0 |
| flavonoid_virtual_bridge | 6 | 0 | 0 | 0 |
| flavonoid_zinc | 6 | 4 | 0 | 0 |

### neohesperidin

- 类别: flavanone_glycoside
- 预期糖家族: neohesperidose

| 实验 | 总路线 | ZINC solved | bridge-closed | strict/trusted |
|------|-------:|------------:|--------------:|---------------:|
| custom_only_virtual_bridge | 5 | 0 | 0 | 0 |
| flavonoid_virtual_bridge | 14 | 0 | 0 | 0 |
| flavonoid_zinc | 13 | 13 | 0 | 0 |

### quercitrin

- 类别: flavonol_glycoside
- 预期糖家族: rhamnoside

| 实验 | 总路线 | ZINC solved | bridge-closed | strict/trusted |
|------|-------:|------------:|--------------:|---------------:|
| custom_only_virtual_bridge | 2 | 0 | 0 | 0 |
| flavonoid_virtual_bridge | 22 | 0 | 0 | 0 |
| flavonoid_zinc | 5 | 0 | 0 | 0 |

### rutin

- 类别: flavonol_glycoside
- 预期糖家族: rutinose

| 实验 | 总路线 | ZINC solved | bridge-closed | strict/trusted |
|------|-------:|------------:|--------------:|---------------:|
| custom_only_virtual_bridge | 2 | 0 | 0 | 0 |
| flavonoid_virtual_bridge | 6 | 0 | 0 | 0 |
| flavonoid_zinc | 8 | 4 | 0 | 0 |

## 聚合统计

### 按实验

| 实验 | 总路线 | ZINC solved | bridge-closed | strict/trusted |
|------|-------:|------------:|--------------:|---------------:|
| custom_only_virtual_bridge | 18 | 0 | 0 | 0 |
| flavonoid_virtual_bridge | 71 | 0 | 11 | 0 |
| flavonoid_zinc | 51 | 38 | 0 | 0 |

### 按糖家族

| 糖家族 | 总路线 | 解决 | 桥接 |
|--------|-------:|-----:|-----:|
| neohesperidose | 59 | 26 | 0 |
| rhamnoside | 29 | 0 | 0 |
| rutinose | 52 | 23 | 11 |
