# 实验框架设计：检测自定义模板和库存的增益

最后更新：2026-07-01

## 核心问题

我们添加的 flavonoid reaction family templates 和分层 stock（strict/trusted/virtual bridge）对原始 USPTO+ZINC 基线是否有显著增益？

## 设计原则

1. **ZINC 是通用库存基线，所有实验都包含 ZINC。** 没有 ZINC 的实验无法回答"增益"问题。
2. **逐层叠加，测量增量。** 每加一个组件，比较前一层的 solved 数。
3. **模板和库存分开消融。** 先测模板增益，再测库存增益。
4. **Panel 为辅，单靶标为主。** 单靶标（hesperidin）做完整消融，panel 做关键实验验证泛化。

## 实验矩阵

### 第一维度：模板增益（库存固定为 ZINC）

| 实验 | Expansion | Stock | 目的 |
|---|---|---|---|
| `A1_uspto_zinc` | USPTO | ZINC | **基线**：标准 USPTO 模板 + 通用库存 |
| `A2_uspto_ringbreaker_zinc` | USPTO + RingBreaker | ZINC | RingBreaker 增益 |
| `A3_uspto_custom_zinc` | USPTO + RingBreaker + Custom | ZINC | **自定义模板增益**（主比较） |

**关键比较：A3 vs A1 → 自定义模板在 ZINC 下的增益**

### 第二维度：库存增益（模板固定为 USPTO+RingBreaker+Custom）

| 实验 | Expansion | Stock | 目的 |
|---|---|---|---|
| `A3_uspto_custom_zinc` | USPTO + RingBreaker + Custom | ZINC | **基线** |
| `B1_plus_strict` | 同上 | ZINC + strict | +strict 增益 |
| `B2_plus_trusted` | 同上 | ZINC + strict + trusted | +trusted 增益 |
| `B3_plus_virtual_bridge` | 同上 | ZINC + strict + trusted + virtual_bridge | +virtual_bridge 增益 |

**关键比较：B3 vs A3 → 分层库存在 ZINC 基础上的增益**

### 第三维度：模板独立性（无 USPTO）

| 实验 | Expansion | Stock | 目的 |
|---|---|---|---|
| `C1_custom_only_zinc` | Custom only | ZINC | 自定义模板独立能力 |
| `C2_custom_only_full_stock` | Custom only | ZINC + strict + trusted + virtual_bridge | 自定义模板 + 全库存 |

**关键比较：C1 vs A1 → 自定义模板能否替代 USPTO**

## 完整实验表

| 编号 | 名称 | Expansion | Stock | 问题 |
|---|---|---|---|---|
| A1 | `uspto_zinc` | USPTO | ZINC | 基线 |
| A2 | `uspto_rb_zinc` | USPTO + RingBreaker | ZINC | RingBreaker 增益？ |
| A3 | `uspto_custom_zinc` | USPTO + RB + Custom | ZINC | Custom 模板增益？ |
| B1 | `uspto_custom_zinc_strict` | 同上 | ZINC + strict | strict 增益？ |
| B2 | `uspto_custom_zinc_trusted` | 同上 | ZINC + strict + trusted | trusted 增益？ |
| B3 | `uspto_custom_zinc_vbridge` | 同上 | ZINC + strict + trusted + VB | VB 增益？ |
| C1 | `custom_only_zinc` | Custom only | ZINC | 模板独立性？ |
| C2 | `custom_only_full_stock` | Custom only | ZINC + strict + trusted + VB | 模板+全库存？ |

## 报告口径

```
增益 = 该实验 solved 数 - 前一层实验 solved 数

模板增益（A3 vs A1）: +X solved → 自定义模板贡献
RingBreaker 增益（A2 vs A1）: +X solved → RingBreaker 贡献
strict 增益（B1 vs A3）: +X solved → strict 库存贡献
trusted 增益（B2 vs B1）: +X solved → trusted 库存贡献
VB 增益（B3 vs B2）: +X solved → virtual bridge 贡献
```

## Panel 实验

Panel 只跑关键实验，不跑全部 8 组：

| Panel 实验 | 对应单靶标实验 | 目的 |
|---|---|---|
| A1 | `uspto_zinc` | ZINC 基线泛化 |
| A3 | `uspto_custom_zinc` | 自定义模板泛化 |
| B3 | `uspto_custom_zinc_vbridge` | 全库存泛化 |

## 路线有效性分类（保持不变）

```
ZINC solved        → 主搜索基准
bridge-closed      → 连通性诊断
strict/trusted     → 高置信子集
unsolved           → 未解决
```

## 文件清单

需要新建/修改的配置：

```
config/ablation_A1_uspto_zinc.yml
config/ablation_A2_uspto_rb_zinc.yml
config/ablation_A3_uspto_custom_zinc.yml
config/ablation_B1_uspto_custom_zinc_strict.yml
config/ablation_B2_uspto_custom_zinc_trusted.yml
config/ablation_B3_uspto_custom_zinc_vbridge.yml
config/ablation_C1_custom_only_zinc.yml
config/ablation_C2_custom_only_full_stock.yml
```

需要修改的脚本：

```
scripts/run_ablation_experiments.py         → 新实验矩阵
scripts/compare_ablation_results.py         → 增益计算逻辑
scripts/run_panel_ablation_experiments.py   → panel 实验对应
scripts/compare_panel_ablation_results.py   → panel 增益报告
```
