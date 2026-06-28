# 橙皮苷专用模板库验证与 tree search 总结报告

## 1. 工作库现状

- 环境：`retro` (Python 3.10, aizynthfinder 4.4.1, RDKit 2023.09.6)
- 模板库：`~/retro_synthesis/templates/flavonoid_templates.hdf5`（清洗后模板）
- 原始 98 个 AI 生成模板 → 自动清洗后 71 个可用，27 个 SMARTS 解析失败
- 已生成的报告：
  - [template_validation_on_hesperidin.csv](~/retro_synthesis/templates/template_validation_on_hesperidin.csv)
  - [flavonoid_templates_failed.csv](~/retro_synthesis/templates/flavonoid_templates_failed.csv)
  - [flavonoid_stock_metadata.csv](~/retro_synthesis/templates/flavonoid_stock_metadata.csv)
  - [flavonoid_stock_inchikeys.txt](~/retro_synthesis/templates/flavonoid_stock_inchikeys.txt)

## 2. 验证结果（清洗后 71 → 增补后 75 个模板）

运行脚本：

```bash
conda run -n retro python /home/ljx/retro_synthesis/src/validate_templates_on_hesperidin.py
```

### 2.1 关键统计

| 指标 | 数值 |
|---|---|
| 模板总数 | 75 |
| 匹配到橙皮苷 | 20 |
| 能裂解出前体 | 20 |
| 裂解产物全在 stock 中 | 0 |

### 2.2 按分类匹配情况

| classification | matched | cracked | total |
|---|---|---|---|
| cyclization | 1 | 1 | 7 |
| deprotection | 14 | 14 | 20 |
| glycosylation | 1 | 1 | 22 |
| oxidation | 1 | 1 | 3 |
| protection | 0 | 0 | 19 |
| reduction | 3 | 3 | 4 |

### 2.3 关键模板状态

- `GLY_009`、`GLY_023`、`DEP_001`、`PROT_004`：产物侧与橙皮苷结构不匹配，未触发匹配。
- `CYC_015`、`GLY_021`：原始 SMARTS 解析失败（在 `flavonoid_templates_failed.csv` 中）。
- 已增补修复版本：
  - `CYC_015_FIXED`：flavanone → chalcone，成功裂解出橙皮苷查尔酮前体。
  - `GLY_021_FIXED`：hesperidin → hesperetin + 游离 rutinose，成功裂解。
  - `DEP_007_FIXED`、`DEP_008_FIXED`、`DEP_013_FIXED`：已加入库中待后续使用。

## 3. 修复内容

1. **stock CSV 修复**：原始 `stock_*.csv` 存在引号未闭合导致 pandas 解析失败的问题。已清洗并重建 `flavonoid_stock_metadata.csv` 与 `flavonoid_stock_inchikeys.txt`。
2. **RDKit / numpy 兼容性补丁**：aizynthfinder 4.4.1 内置 `DataStructs.ConvertToNumpyArray` 与 numpy 2.x 不兼容，已修补 `aizynthfinder/chem/mol.py` 中的 fingerprint 生成逻辑为 `np.array(bitvect, dtype=np.float32)`。
3. **模板增补**：通过 `augment_templates.py` 将 4 个手工修复模板加入 `flavonoid_templates.hdf5` / `.csv.gz`，总数从 71 提升到 75。
4. **无模型自定义扩展策略**：为避免专用模板数量（75）与 USPTO 策略网络输出维度不匹配，编写了 `custom_expansion.UniformTemplateExpansion`，对匹配的模板均匀分配先验概率，无需训练专用策略网络即可调用。

## 4. Tree search 配置与结果

### 4.1 配置文件

- 组合策略：`~/retro_synthesis/config/flavonoid.yml`
  - expansion：`uspto`（USPTO 策略网络 + 42554 模板） + `flavonoid`（UniformTemplateExpansion，仅 glycosylation/cyclization 两类）
  - stock：`zinc` + `flavonoid`
  - filter：`uspto_filter_model.onnx`
  - search：MCTS，500 iterations，600 s，max_transforms=10

### 4.2 组合搜索（USPTO + 专用模板）

```bash
conda run -n retro python -m aizynthfinder.interfaces.aizynthcli \
  --smiles "COc1ccc(C2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)[C@H](O)[C@@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O" \
  --config config/flavonoid.yml --policy uspto flavonoid --stocks zinc flavonoid \
  --output outputs/hesperidin_combined_curated.json --log_to_file
```

结果：

- 找到 6 条 route，is_solved: True
- Top route：2 步，先甲基化再打开 flavanone 环得到查尔酮（已在 stock 中）
- 其他 route：6-7 步，走 USPTO 的 flavanone C-C 键断裂 → isovanillin + phloracetophenone 糖苷路线，随后进行保护/去保护/糖基切断。

### 4.3 USPTO 单独搜索

```bash
conda run -n retro python -m aizynthfinder.interfaces.aizynthcli \
  --smiles "..." --config config/flavonoid.yml --policy uspto --stocks zinc flavonoid \
  --output outputs/hesperidin_uspto_only.json --log_to_file
```

结果：

- 找到 5 条 solved route，4-6 步，全部前体在 stock 中。
- 路线同样以 flavanone C-C 断裂 + 糖上保护/去保护/氟代切断为主，化学可解释性一般。

### 4.4 专用模板单独搜索

```bash
conda run -n retro python -m aizynthfinder.interfaces.aizynthcli \
  --smiles "..." --config config/flavonoid.yml --policy flavonoid --stocks flavonoid \
  --output outputs/hesperidin_custom_only.json --log_to_file
```

结果：

- 仅 4 个 node，2 条 route，未 solved。
- `GLY_021_FIXED` 成功将橙皮苷裂解为 hesperetin（在 stock）+ rutinose 片段（不在 stock，因 RDKit 反应输出为非常规糖环表示）。

## 5. 与 USPTO 通用模型的对比

| 维度 | USPTO 通用模型 | 专用模板库（当前） |
|---|---|---|
| 模板数量 | 42,554 | 75 |
| 策略网络 | 已训练 | 暂无（uniform prior） |
| 对橙皮苷的 top route | 5-6 步，走通用官能团转换 | 可 1-2 步直接切到查尔酮/hesperetin |
| 化学合理性 | 一般（出现保护基乱用、氟代切糖） | 需大量清洗，方向/芳香性错误多 |
| 运行速度 | 较快 | 较慢（uniform 扩展产生较多候选） |
| 可解释性 | 低 | 可针对黄酮/糖苷反应人工设计 |

结论：专用模板库在概念上更适合橙皮苷，但目前 27 个失败模板和过于宽泛的保护/去保护模板导致搜索噪声大；USPTO 虽能找到 route，但多是数据驱动的非常规断键。两者结合时，USPTO 的评分机制仍占主导，专用模板的优势未充分发挥。

## 6. 进一步优化方向

1. **模板精修**
   - 逐个修复 27 个失败模板中的关键反应（`CYC_015`、`GLY_021`、`DEP_007/008/013` 等）。
   - 删除或限定过于通用的保护/去保护模板，避免在搜索中无限添加 TBS/Bn/PMB。
   - 为 hesperidin 建立专用 glycosylation 模板：flavonoid-O-rutinoside → flavonoid-OH + 活化 rutinose donor。

2. **专用策略网络**
   - 基于现有 75 个模板训练一个小的分类器/策略网络（输出 75 维），让 tree search 学会优先选择黄酮/糖苷相关模板。
   - 数据增强：从文献或 USPTO 中提取更多 hesperidin 相关反应实例，扩展训练集。

3. **Stock 与评分优化**
   - 扩充专用 stock：加入 hesperetin、naringenin chalcone、活化糖基供体、常见保护基试剂。
   - 引入 reaction class scorer / cost scorer，惩罚保护基循环和不自然断键。
   - 使用 `break_bonds` / `freeze_bonds` 约束，强制优先切断糖苷键和 flavanone C-ring。

4. **多目标/约束搜索**
   - 设置最长步数、保护基计数、分子复杂度等约束。
   - 结合 Retro* 或 MCTS 的 custom reward，鼓励 route 终止于已知的黄酮前体。

5. **验证闭环**
   - 对每条 top route 用 RDKit 正向反应验证，确保 SMARTS 方向、原子映射、立体化学正确。
   - 引入人工化学审查，淘汰 USPTO 中产生的保护基/卤代等不合理步骤。

## 7. 生成的关键文件

| 文件 | 说明 |
|---|---|
| `templates/template_validation_on_hesperidin.csv` | 75 个模板对橙皮苷的匹配/裂解报告 |
| `templates/flavonoid_templates.hdf5` | 增补后的 aizynthfinder 模板库 |
| `templates/flavonoid_templates.csv.gz` | 带元数据的完整模板表 |
| `templates/flavonoid_stock_inchikeys.txt` | 专用 stock InChIKey 列表 |
| `config/flavonoid.yml` | USPTO + 专用模板的 tree search 配置 |
| `scripts/custom_expansion.py` | 无模型 uniform 扩展策略 |
| `outputs/hesperidin_combined_curated.json` | 组合搜索 top routes |
| `outputs/hesperidin_uspto_only.json` | USPTO 单独搜索 routes |
| `outputs/hesperidin_custom_only.json` | 专用模板单独搜索 routes |
