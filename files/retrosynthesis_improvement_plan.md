# Aizynthfinder 针对黄酮苷类化合物的逆合成分析能力提升计划

## 背景与目标

Aizynthfinder 基于 USPTO 通用数据库构建，对于黄酮苷类化合物（以橙皮苷 hesperidin 为代表）的逆合成分析准确度不足，尤其对 **O-糖苷键**的切割策略识别不明确。本计划旨在通过针对性扩充 reaction templates、stock、以及模型微调/参数优化，显著提升逆合成路径可行性。

---

## 第一阶段：问题诊断与基准建立（第1-2周）

### 1.1 Hesperidin 结构解剖

橙皮苷 = 橙皮素(hesperetin) + 芸香糖(rutinose, α-L-鼠李糖-(1→6)-β-D-葡萄糖)

需要关注的关键断键位点：
- **O-糖苷键**：芸香糖 C1 与橙皮素 7-OH 之间的连接（最关键的瓶颈）
- **糖内 C-O-C 键**：鼠李糖与葡萄糖之间的 (1→6) 连接
- **C-C 键**：黄酮骨架内的可能逆合成断键（C-ring 开环等）

### 1.2 基线性能评估

在当前 USPTO 模板集上运行以下测试：

| 测试化合物 | 结构特征 | 预期断键 |
|---|---|---|
| Hesperidin (橙皮苷) | flavanone-7-O-rutinoside | O-糖苷键 |
| Rutin (芦丁) | flavonol-3-O-rutinoside | O-糖苷键 |
| Naringin (柚皮苷) | flavanone-7-O-neohesperidoside | O-糖苷键 |
| Quercitrin (槲皮苷) | flavonol-3-O-rhamnoside | O-糖苷键 |
| Vitexin (牡荆素) | flavone-8-C-glucoside | C-糖苷键（对照） |

记录每条逆合成路径的：
- 每步使用的 reaction template
- 路径总步数
- 每步的可行性评分（policy network 分数）
- 最终能否到达 stock 中的起始原料

### 1.3 失败原因分析

典型失败模式预判：
1. **模板缺失**：USPTO 中 O-糖苷化反应的覆盖不足，尤其是黄酮类底物
2. **Stock 缺失**：橙皮素、芸香糖、鼠李糖等黄酮/糖类构建模块不在 stock 中
3. **策略偏差**：模型可能选择切断 C-ring 或 A/B-ring 等非最优路径，而非最直接的 O-糖苷键切割
4. **立体化学**：糖苷键的 α/β 构型未被充分区分

---

## 第二阶段：Reaction Template 扩充（第3-6周）

### 2.1 数据来源

| 来源 | 内容 | 优先级 |
|---|---|---|
| **USPTO 原始数据** | 重新筛选含糖苷化反应的条目，提取更细粒度的模板 | 高 |
| **Reaxys / SciFinder** | 黄酮苷合成相关文献反应 | 高 |
| **Brenda / KEGG** | 生物催化（糖基转移酶）反应 | 中 |
| **文献手编** | 经典糖苷化方法学（Koenigs-Knorr, Schmidt, 三氯乙酰亚胺酯法等） | 中 |

### 2.2 模板提取策略

针对 O-糖苷键切割，提取以下类别的 reaction templates：

**A. 糖苷化反应（逆向 = 糖苷键断裂）**
- 羟基糖供体 + 糖受体 → 糖苷（经典糖苷化）
- 酶催化糖基转移反应（UDP-糖供体 → 糖苷）
- 一锅法区域选择性糖苷化

**B. 保护/去保护策略**
- 糖羟基的选择性保护（乙酰化、苄基化、硅醚化）
- 黄酮羟基的保护/去保护

**C. 黄酮骨架构建（如需逆合成到更简单的前体）**
- Baker-Venkataraman 重排
- Allan-Robinson 反应
- 查尔酮环化

### 2.3 模板格式化

遵循 aizynthfinder 的 SMARTS-based template 格式：
```
{
  "reaction_smarts": "[#6:1]-[OX2:2]-[c:3]>>[#6:1]-[OX2:2].[c:3]-[OX2:4]-[CH]...", 
  "reference": "文献来源",
  "category": "glycosylation",
  "confidence": 0.85,
  "metadata": {
    "reaction_type": "O-glycosylation",
    "substrate_class": "flavonoid",
    "typical_yield": "60-85%"
  }
}
```

关键：模板需要包含适当的反应中心上下文（context radius），以区分：
- O-糖苷键 vs C-糖苷键
- 不同位置的羟基糖苷化（3-O, 7-O, 4'-O）
- α vs β 糖苷键

### 2.4 验证新模板

用 RDKit 检查每条新模板的：
- SMARTS 语法正确性
- 产物可映射性（反应后原子映射一致性）
- 与现有模板无冲突/重复

---

## 第三阶段：Stock 库扩充（第3-5周，与第二阶段并行）

### 3.1 核心糖类构建模块

| 分子 | SMILES（示意） | 用途 |
|---|---|---|
| β-D-Glucose (带适当保护基) | - | 葡萄糖供体 |
| α-L-Rhamnose | - | 鼠李糖供体 |
| Rutinose (芸香糖) | - | 整体供体（若可商业化） |
| UDP-glucose | - | 生物催化供体 |
| UDP-rhamnose | - | 生物催化供体 |

### 3.2 黄酮类起始原料

| 分子 | 说明 |
|---|---|
| Hesperetin (橙皮素) | hesperidin 的糖苷配基 |
| Naringenin (柚皮素) | flavanone 骨架 |
| Quercetin (槲皮素) | flavonol 骨架 |
| 2,4,6-三羟基苯乙酮 | Baker-Venkataraman 前体 |
| 对羟基苯甲醛衍生物 | 查尔酮前体 |

### 3.3 Stock 来源

- Sigma-Aldrich / TCI 等商业可得化合物（优先）
- PubChem CID 查询确认可得性
- 对于非商业化化合物，标记为"文献合成可得"并添加合成路线

---

## 第四阶段：模型优化（第6-10周）

### 4.1 策略选择

有两种可行的优化路径，可并行或按优先级推进：

#### 路径 A：模板库增强 + 参数调优（推荐先行）

无需重新训练神经网络，通过以下方式优化：

1. **调整模板库组成**
   - 将新提取的黄酮/糖苷化模板整合进 `templates.json`
   - 调整 template frequency 上的置信度权重

2. **Stock 注释优化**
   - 将扩充的 stock 以正确的格式写入 `stock.json`
   - 添加 molecular fingerprint 索引以加速子结构匹配

3. **搜索参数调优**
   - `max_transforms`：增加每步考虑的模板数量（默认50 → 100-150）
   - `min_history`：调整剪枝阈值
   - `cutoff_cumulative` / `cutoff_number`：放宽早停条件
   - 搜索算法选择：UCB1 vs PUCT 参数（C值调整）

4. **过滤器网络参数**
   - 若使用 filter policy，调整 filter threshold
   - 对黄酮类底物适当放宽过滤条件

#### 路径 B：神经网络微调（更高成本，更显著提升）

1. **训练数据准备**
   - 从 Reaxys/USPTO 筛选 flavonoid glycoside 相关反应（约2000-5000条）
   - 化学反应 SMILES → template 提取
   - 划分 train/val/test（80/10/10）

2. **Expansion Policy Network 微调**
   - 以现有权重为初始化（迁移学习）
   - 小学习率微调（lr = 1e-4 ~ 1e-5）
   - Epoch: 20-50，早停监控 validation loss
   - 关键指标：glycosylation template 的召回率

3. **Filter Policy Network 微调（如适用）**
   - 正样本：可行的黄酮苷逆合成步骤
   - 负样本：不可行的断键策略
   - 关注 precision@k（前 k 个推荐中可行的比例）

### 4.2 技术实现细节

```
aizynthfinder 版本: 需确认 (建议 ≥ 3.x)
数据格式: JSON-based templates, HDF5 model weights
依赖: RDKit, TensorFlow/PyTorch (取决于版本), numpy
```

---

## 第五阶段：验证与评估（第10-12周）

### 5.1 定量评估

| 指标 | 定义 | 基线目标 |
|---|---|---|
| **Route Finding Rate** | 成功找到完整逆合成路径的比例 | 提升至 >80% |
| **Feasibility Score** | 路径中每步可行性的平均分 | 显著提升 |
| **Template Hit Rate** | 遇到 O-糖苷键时正确选择 glycosylation template 的比例 | 从基线提升 |
| **Average Steps** | 完成路径的平均步数 | 合理范围 (3-8步) |
| **Stock Reachability** | 最终前体落在 stock 中的比例 | 提升至 >90% |

### 5.2 定性评估（专家评审）

- 请合成化学专家评审 top-3 推荐路径的可行性
- 重点关注：
  - 化学选择性是否合理（哪些位点先反应）
  - 保护/去保护策略是否必要且正确
  - 立体化学控制是否被考虑

### 5.3 消融实验

| 实验 | 变量 | 目的 |
|---|---|---|
| Baseline | 原始 USPTO templates + stock | 对照 |
| +Templates | 仅扩充 templates | 模板贡献 |
| +Stock | 仅扩充 stock | stock 贡献 |
| +Templates +Stock | 同时扩充 | 交互效应 |
| +Fine-tuning | 完整流程 | 整体提升 |

---

## 第六阶段：扩展与泛化（第12周+）

- 将方法扩展到更广泛的 flavonoid glycoside 子类
- 考虑 C-糖苷键（如 vitexin）的策略
- 整理为可复用的 pipeline，为后续其他天然产物类做准备

---

## 时间线概览

```
Week 1-2    ████████  问题诊断 & 基线建立
Week 3-6    ████████████████  Template 扩充
Week 3-5    ████████████  Stock 扩充
Week 6-10   ████████████████████  模型优化 (A: 参数调优 / B: 微调)
Week 10-12  ████████  验证 & 评估
Week 12+    ████  扩展 & 泛化
```

---

## 关键风险与缓解

| 风险 | 影响 | 缓解措施 |
|---|---|---|
| 手动模板质量参差 | 产生不可行的路径 | RDKit 严格验证 + 专家审查 |
| 训练数据不足（flavonoid 反应有限） | 微调效果受限 | 优先走路径 A (参数调优)；路径 B 作为补充 |
| 模板冲突（新旧模板矛盾） | 降低整体搜索效率 | 版本控制模板库，逐步集成并 A/B 测试 |
| aizynthfinder 版本差异 | 代码接口变化 | 早期确认版本，固定依赖 |

---

## 下一步行动

1. 确认当前使用的 aizynthfinder 版本和运行环境
2. 对 hesperidin 运行一次完整的基线分析，记录失败模式
3. 开始 USPTO 数据中 glycosylation 反应的筛选和模板提取
4. 搭建消融实验的评估框架
