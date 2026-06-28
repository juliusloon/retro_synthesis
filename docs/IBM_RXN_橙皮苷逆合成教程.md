# IBM RXN for Chemistry 橙皮苷逆合成完整教程

## 前置条件

- Python 3.8+
- IBM RXN API Key（已从 IBM RXN 网站获取）
- 网络连接（可访问 IBM RXN 云服务）

---

## 第一步：环境准备

### 1.1 创建虚拟环境

```bash
# 创建虚拟环境
python -m venv rxn_env

# 激活（Linux/Mac）
source rxn_env/bin/activate

# 激活（Windows）
rxn_env\Scripts\activate
```

### 1.2 安装依赖包

```bash
pip install rxn4chemistry rdkit requests
```

验证安装：

```python
import rxn4chemistry
print(rxn4chemistry.__version__)
```

---

## 第二步：获取 API Key 并初始化

### 2.1 获取项目 ID

IBM RXN 使用项目（Project）来组织任务。首次使用需要创建项目。

```python
from rxn4chemistry import RXN4ChemistryWrapper

# 替换为你的实际 API Key
API_KEY = "your-ibm-rxn-api-key-here"

# 初始化客户端
rxn = RXN4ChemistryWrapper(api_key=API_KEY)

# 创建新项目（如果已有项目 ID 可跳过）
response = rxn.create_project(name="Hesperidin_Retrosynthesis")
print(response)

# 记录项目 ID
PROJECT_ID = response['response']['payload']['id']
print(f"项目 ID: {PROJECT_ID}")
```

### 2.2 查看已有项目（可选）

```python
# 列出所有项目
projects = rxn.list_all_projects()
# 注意：项目列表在 payload['content'] 中，不是 payload['projects']
for project in projects['payload']['content']:
    print(f"项目名: {project['name']}, ID: {project['id']}")
```

---

## 第三步：准备橙皮苷分子

### 3.1 橙皮苷 SMILES

橙皮苷（Hesperidin）的 SMILES 表示：

```
COc1ccc(C2CC(=O)c3c(cc(OC4OC(COC5OC(C)C(O)C(O)C5O)C(O)C(O)C4O)c(OC)c3O)O2)cc1O
```

### 3.2 用 RDKit 验证结构

```python
from rdkit import Chem
from rdkit.Chem import Draw

hesperidin_smiles = "COc1ccc(C2CC(=O)c3c(cc(OC4OC(COC5OC(C)C(O)C(O)C5O)C(O)C(O)C4O)c(OC)c3O)O2)cc1O"

mol = Chem.MolFromSmiles(hesperidin_smiles)
if mol is None:
    print("SMILES 解析失败！")
else:
    print("SMILES 解析成功")
    print(f"分子式: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
    print(f"分子量: {Chem.rdMolDescriptors.CalcExactMolWt(mol):.2f}")
    
    # 保存分子图片
    img = Draw.MolToImage(mol, size=(400, 400))
    img.save("hesperidin_structure.png")
    print("分子结构图已保存为 hesperidin_structure.png")
```

---

## 第四步：执行逆合成预测

### 4.1 提交逆合成任务

```python
from rxn4chemistry import RXN4ChemistryWrapper
import time

API_KEY = "your-ibm-rxn-api-key-here"
PROJECT_ID = "your-project-id-here"

rxn = RXN4ChemistryWrapper(api_key=API_KEY)
rxn.set_project(PROJECT_ID)

hesperidin_smiles = "COc1ccc(C2CC(=O)c3c(cc(OC4OC(COC5OC(C)C(O)C(O)C5O)C(O)C(O)C4O)c(OC)c3O)O2)cc1O"

# 提交逆合成预测请求
print("正在提交逆合成预测请求...")
response = rxn.predict_automatic_retrosynthesis(
    product=hesperidin_smiles,
    max_steps=5,           # 最大步数
    nbeams=10,             # 搜索束宽
    pruning_steps=2,       # 剪枝步数
    ai_model="2020-07-01"  # 模型版本
)

# 获取预测任务 ID
prediction_id = response['response']['payload']['id']
print(f"预测任务 ID: {prediction_id}")

# 保存任务 ID 供后续查询
with open("prediction_id.txt", "w") as f:
    f.write(prediction_id)
```

### 4.2 关键参数说明

| 参数 | 说明 | 建议值 |
|------|------|--------|
| `product` | 目标产物 SMILES | 橙皮苷 SMILES |
| `max_steps` | 最大逆推步数 | 3–5（橙皮苷结构复杂，建议 5） |
| `nbeams` | 搜索束宽 | 10–20（越大搜索越充分） |
| `pruning_steps` | 剪枝步数 | 2–3（控制搜索空间） |
| `ai_model` | AI 模型版本 | "2020-07-01" 或最新版 |

---

## 第五步：轮询等待结果

逆合成是异步计算，需要轮询等待完成。

```python
import time

print("等待逆合成计算完成...")
print("通常需要 30 秒到 5 分钟，取决于分子复杂度")

max_wait = 600  # 最大等待 10 分钟
waited = 0
interval = 10   # 每 10 秒查询一次

while waited < max_wait:
    result = rxn.get_predict_automatic_retrosynthesis_results(prediction_id)
    status = result['response']['payload']['status']
    
    print(f"[{waited}s] 状态: {status}")
    
    if status == 'SUCCESS':
        print("计算完成！")
        break
    elif status in ['FAILED', 'ERROR']:
        print(f"计算失败: {result}")
        break
    
    time.sleep(interval)
    waited += interval
else:
    print("等待超时，请稍后手动查询结果")

# 保存完整结果
import json
with open("retrosynthesis_result.json", "w", encoding="utf-8") as f:
    json.dump(result, f, ensure_ascii=False, indent=2)
print("结果已保存到 retrosynthesis_result.json")
```

---

## 第六步：解析和展示结果

### 6.1 提取合成路径

```python
import json

with open("retrosynthesis_result.json", "r", encoding="utf-8") as f:
    result = json.load(f)

payload = result['response']['payload']

# 检查是否有结果
if 'retrosynthetic_paths' not in payload or not payload['retrosynthetic_paths']:
    print("未找到逆合成路径，可能原因：")
    print("- 分子过于复杂，超出模型能力范围")
    print("- 模型未覆盖相关反应类型（如糖基化）")
    print("建议尝试更简单的分子或调整参数")
else:
    paths = payload['retrosynthetic_paths']
    print(f"共找到 {len(paths)} 条逆合成路径\n")
    
    for i, path in enumerate(paths[:5], 1):  # 只展示前 5 条
        print(f"=== 路径 {i} ===")
        print(f"置信度: {path.get('confidence', 'N/A')}")
        print(f"步数: {len(path.get('reactions', []))}")
        
        # 遍历每一步反应
        for step_num, reaction in enumerate(path.get('reactions', []), 1):
            print(f"\n  步骤 {step_num}:")
            print(f"    反应 SMILES: {reaction.get('smiles', 'N/A')}")
            print(f"    反应名称: {reaction.get('name', 'N/A')}")
            print(f"    置信度: {reaction.get('confidence', 'N/A')}")
            
            # 反应物
            reactants = reaction.get('reactants', [])
            print(f"    反应物: {', '.join(reactants)}")
            
            # 产物
            products = reaction.get('products', [])
            print(f"    产物: {', '.join(products)}")
        
        print("\n" + "-" * 50)
```

### 6.2 解析 IBM RXN 的树状路径

IBM RXN 返回的是**嵌套树结构**（每个节点有 `children`），不是扁平列表。需要先递归展开成线性路线：

```python
import json
from rdkit import Chem
from rdkit.Chem import Draw

with open("result.json", "r", encoding="utf-8") as f:
    result = json.load(f)

# 路径在顶层 key 'retrosynthetic_paths' 中
path_trees = result.get('retrosynthetic_paths', [])
print(f"共 {len(path_trees)} 棵路径树")

def extract_routes(node, current_path=None):
    """递归提取所有从根到叶子的线性路径
    
    每个节点包含:
      - smiles: 当前分子的 SMILES
      - confidence: 该步的置信度
      - isCommercial: 是否为商业可得原料
      - children: 下一层反应物列表
    """
    if current_path is None:
        current_path = []
    
    current_path.append({
        'smiles': node.get('smiles', ''),
        'confidence': node.get('confidence', 0),
        'isCommercial': node.get('isCommercial', False),
        'rclass': node.get('rclass', 'Undefined')  # 反应类型
    })
    
    children = node.get('children', [])
    if not children:
        # 到达叶子节点，返回完整路径
        yield list(current_path)
    else:
        for child in children:
            yield from extract_routes(child, list(current_path))

# 展开所有路径
all_routes = []
for tree in path_trees:
    for route in extract_routes(tree):
        all_routes.append(route)

print(f"共解析出 {len(all_routes)} 条完整合成路线\n")

# 展示前3条路线
for i, route in enumerate(all_routes[:3], 1):
    print(f"=== 路线 {i} ({len(route)} 步) ===")
    for j, step in enumerate(route):
        marker = " ✓" if step['isCommercial'] else ""
        print(f"  步骤 {j}: {step['smiles']}{marker}")
    print()
```

### 6.3 可视化单步反应

IBM RXN 的树结构中，**每一步逆合成本质是**：一个产物（父节点）被拆分成若干反应物（子节点）。可视化时需要把这一步表示为 "产物 >> 反应物" 的反应式。

```python
def visualize_retro_step(parent_smiles, children_nodes, step_num, output_dir="."):
    """可视化逆合成单步: 产物 >> 反应物1 + 反应物2 + ...
    
    Args:
        parent_smiles: 产物 SMILES（当前节点）
        children_nodes: 子节点列表（反应物）
        step_num: 步骤编号（用于文件名）
        output_dir: 输出目录
    """
    try:
        # 构造反应 SMARTS: 产物 >> 反应物
        reactants = ".".join([child['smiles'] for child in children_nodes])
        reaction_smarts = f"{parent_smiles}>>{reactants}"
        
        # 尝试解析反应
        rxn = Chem.rdChemReactions.ReactionFromSmarts(reaction_smarts, useSmiles=True)
        if rxn is None:
            # 如果 SMARTS 解析失败，分别画分子拼图
            print(f"步骤 {step_num}: 使用分屏可视化")
            molecules = [Chem.MolFromSmiles(parent_smiles)] + \
                       [Chem.MolFromSmiles(c['smiles']) for c in children_nodes]
            labels = ["产物"] + [f"反应物{i+1}" for i in range(len(children_nodes))]
            img = Draw.MolsToGridImage(molecules, molsPerRow=3, subImgSize=(300, 300),
                                       legends=labels)
            filename = f"{output_dir}/step_{step_num}_grid.png"
        else:
            img = Draw.ReactionToImage(rxn, subImgSize=(400, 300))
            filename = f"{output_dir}/step_{step_num}_reaction.png"
        
        img.save(filename)
        print(f"步骤 {step_num} 已保存: {filename}")
        return filename
        
    except Exception as e:
        print(f"步骤 {step_num} 可视化失败: {e}")
        return None

# 可视化第一条路线的每一步
if all_routes:
    route = all_routes[0]
    print(f"\n可视化路线 1（共 {len(route)-1} 步逆合成）：\n")
    
    # 递归地可视化每一步：从根往下，每个有 children 的节点就是一步
    def visualize_route_recursive(node, step=1):
        children = node.get('children', [])
        if children:
            # 这是一步逆合成：当前节点 = 产物，children = 反应物
            visualize_retro_step(node['smiles'], children, step)
            # 递归处理每个子节点的下一步
            for child in children:
                next_step = visualize_route_recursive(child, step + 1)
                step = max(step, next_step)
        return step
    
    # 从第一棵树的根开始可视化
    visualize_route_recursive(path_trees[0])
```

### 6.4 提取起始原料

起始原料是路径树中的**叶子节点**（没有 children），通常标记为 `isCommercial=true`：

```python
def extract_starting_materials(route):
    """从一条线性路径中提取起始原料
    
    返回: list of dict，每个包含 smiles, confidence, isCommercial
    """
    # 路线的最后一个节点是叶子 = 起始原料
    # 但实际上逆合成树的叶子就是商业可得的起始物
    materials = []
    for step in route:
        if step['isCommercial']:
            materials.append(step)
    return materials

# 分析每条路线的起始原料
print("\n=== 各路线起始原料分析 ===\n")
for i, route in enumerate(all_routes[:5], 1):
    materials = extract_starting_materials(route)
    print(f"路线 {i}: {len(materials)} 个起始原料")
    
    for mat in materials:
        mol = Chem.MolFromSmiles(mat['smiles'])
        mw = Chem.rdMolDescriptors.CalcExactMolWt(mol) if mol else 0
        print(f"  - {mat['smiles']}")
        print(f"    分子量: {mw:.1f}, 置信度: {mat['confidence']}")
    print()

# 统计所有路线的独特起始原料
all_materials = set()
for route in all_routes:
    for step in route:
        if step['isCommercial']:
            all_materials.add(step['smiles'])

print(f"所有路线共涉及 {len(all_materials)} 种不同的起始原料")
for smiles in sorted(all_materials)[:10]:  # 展示前10个
    print(f"  - {smiles}")
```

---

## 第七步：高级用法

### 7.1 指定特定反应类型（过滤）

IBM RXN 支持对特定化学空间进行限制：

```python
# 请求中可添加约束（需查阅最新 API 文档）
# 例如：限制保护基反应的使用
response = rxn.predict_automatic_retrosynthesis(
    product=hesperidin_smiles,
    max_steps=5,
    nbeams=10,
    # 可添加其他高级参数
)
```

### 7.2 批量预测多个分子

```python
molecules = [
    ("Hesperidin", "COc1ccc(C2CC(=O)c3c(cc(OC4OC(COC5OC(C)C(O)C(O)C5O)C(O)C(O)C4O)c(OC)c3O)O2)cc1O"),
    ("Hesperetin", "COc1ccc(C2CC(=O)c3c(cc(O)c(OC)c3O)O2)cc1O"),  # 橙皮素（苷元）
    ("Glucose", "C(C1C(C(C(C(O1)O)O)O)O)O"),  # 葡萄糖（简单验证）
]

results = {}
for name, smiles in molecules:
    print(f"\n处理: {name}")
    response = rxn.predict_automatic_retrosynthesis(product=smiles, max_steps=3, nbeams=5)
    pred_id = response['response']['payload']['id']
    results[name] = pred_id
    print(f"任务 ID: {pred_id}")
    time.sleep(2)  # 避免请求过快

# 保存所有任务 ID
with open("batch_predictions.json", "w") as f:
    json.dump(results, f, indent=2)
```

### 7.3 使用正合成验证逆合成路线

```python
# 选择一条逆合成路径，用正合成验证其可行性
def validate_route_forward(path):
    """对逆合成路径进行正合成验证"""
    reactions = path.get('reactions', [])
    
    # 从后往前（逆合成的顺序），正向验证每一步
    for i, reaction in enumerate(reversed(reactions), 1):
        reactants = ".".join(reaction.get('reactants', []))
        print(f"验证步骤 {i}: {reactants} >> {reaction.get('products', [''])[0]}")
        
        # 调用正合成预测
        fwd_response = rxn.predict_reaction(reactants)
        print(f"正合成预测结果: {fwd_response}")
        time.sleep(1)

# 验证第一条路径
if paths:
    validate_route_forward(paths[0])
```

---

## 第八步：常见问题与解决

### Q1: API 返回 "Invalid API Key"

```python
# 检查 API Key 是否正确
print(f"API Key 长度: {len(API_KEY)}")
print(f"API Key 前缀: {API_KEY[:8]}...")

# 重新初始化
rxn = RXN4ChemistryWrapper(api_key=API_KEY)
```

### Q2: 橙皮苷计算超时或失败

橙皮苷是大分子（分子量 ~610），糖基化结构复杂，IBM RXN 可能难以处理。

**解决方案：**
1. 先测试更简单的分子（如葡萄糖、苯甲酸）验证 API 可用
2. 尝试拆分分析：分别预测橙皮素（苷元）和糖基部分的合成
3. 调整参数：增大 nbeams，减小 max_steps

```python
# 测试简单分子
test_smiles = "c1ccccc1C(=O)O"  # 苯甲酸
response = rxn.predict_automatic_retrosynthesis(product=test_smiles, max_steps=3)
print("简单分子测试:", response)
```

### Q3: 结果为空或无路径

可能原因：
- 分子包含模型未训练过的结构（如特定糖苷键）
- 参数过于严格（pruning_steps 过大）

**解决方案：**
```python
# 放宽搜索参数
response = rxn.predict_automatic_retrosynthesis(
    product=hesperidin_smiles,
    max_steps=5,
    nbeams=20,        # 增大搜索空间
    pruning_steps=1   # 减少剪枝
)
```

### Q4: 如何获取反应条件信息

IBM RXN 的逆合成 API 主要返回反应物和产物，**不直接提供反应条件**。但你可以：
1. 通过反应名称在 Reaxys 或 SciFinder 中查找具体条件
2. 使用 IBM RXN 的正合成 API 获取更详细的预测结果

---

## 第九步：完整脚本

以下是一个可直接运行的完整脚本：

```python
#!/usr/bin/env python3
"""
IBM RXN 橙皮苷逆合成分析完整脚本
运行前请替换 API_KEY 和 PROJECT_ID
"""

import time
import json
from rxn4chemistry import RXN4ChemistryWrapper
from rdkit import Chem
from rdkit.Chem import Draw

# ==================== 配置 ====================
API_KEY = "your-ibm-rxn-api-key-here"
PROJECT_ID = None  # 设为 None 则自动创建新项目

HESPERIDIN_SMILES = "COc1ccc(C2CC(=O)c3c(cc(OC4OC(COC5OC(C)C(O)C(O)C5O)C(O)C(O)C4O)c(OC)c3O)O2)cc1O"

# ==================== 初始化 ====================
rxn = RXN4ChemistryWrapper(api_key=API_KEY)

if PROJECT_ID is None:
    response = rxn.create_project(name="Hesperidin_Analysis")
    PROJECT_ID = response['response']['payload']['id']
    print(f"创建新项目: {PROJECT_ID}")
else:
    print(f"使用已有项目: {PROJECT_ID}")

rxn.set_project(PROJECT_ID)

# ==================== 验证分子 ====================
mol = Chem.MolFromSmiles(HESPERIDIN_SMILES)
print(f"\n目标分子: 橙皮苷 (Hesperidin)")
print(f"分子式: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
print(f"分子量: {Chem.rdMolDescriptors.CalcExactMolWt(mol):.2f}")

# 保存结构图
img = Draw.MolToImage(mol, size=(400, 400))
img.save("hesperidin.png")
print("结构图已保存: hesperidin.png")

# ==================== 提交逆合成请求 ====================
print("\n提交逆合成预测...")
response = rxn.predict_automatic_retrosynthesis(
    product=HESPERIDIN_SMILES,
    max_steps=5,
    nbeams=10,
    pruning_steps=2
)

prediction_id = response['response']['payload']['id']
print(f"任务 ID: {prediction_id}")

# ==================== 等待结果 ====================
print("\n等待计算完成...")
max_wait = 600
waited = 0
interval = 10

while waited < max_wait:
    result = rxn.get_predict_automatic_retrosynthesis_results(prediction_id)
    status = result['response']['payload']['status']
    print(f"  [{waited}s] 状态: {status}")
    
    if status == 'SUCCESS':
        break
    elif status in ['FAILED', 'ERROR']:
        print(f"计算失败!")
        exit(1)
    
    time.sleep(interval)
    waited += interval

# ==================== 保存结果 ====================
with open("result.json", "w", encoding="utf-8") as f:
    json.dump(result, f, ensure_ascii=False, indent=2)

# ==================== 解析和展示 ====================
payload = result['response']['payload']
paths = payload.get('retrosynthetic_paths', [])

if not paths:
    print("\n未找到逆合成路径")
    print("建议: 尝试更简单的分子，或调整 max_steps/nbeams 参数")
else:
    print(f"\n✓ 找到 {len(paths)} 条路径")
    
    for i, path in enumerate(paths[:3], 1):
        print(f"\n{'='*60}")
        print(f"路径 {i} / 置信度: {path.get('confidence', 'N/A')}")
        print(f"{'='*60}")
        
        for step_num, reaction in enumerate(path.get('reactions', []), 1):
            print(f"\n步骤 {step_num}:")
            print(f"  反应: {reaction.get('smiles', 'N/A')}")
            print(f"  名称: {reaction.get('name', 'N/A')}")
            print(f"  反应物: {reaction.get('reactants', [])}")
            print(f"  产物: {reaction.get('products', [])}")

print("\n✓ 分析完成，详细结果保存在 result.json")
```

---

## 参考资源

- IBM RXN 官方文档: https://rxn.res.ibm.com/
- rxn4chemistry Python 包: https://github.com/rxn4chemistry/rxn4chemistry
- IBM RXN API 文档: https://rxn.res.ibm.com/rxn/apidoc/
- 橙皮苷 PubChem 页面: https://pubchem.ncbi.nlm.nih.gov/compound/Hesperidin

---

## 下一步建议

1. **跑通本教程**：先验证 API 和基本流程是否正常
2. **对比分析**：将 IBM RXN 的结果与你的 AiZynthFinder 结果对比
3. **提取模板**：从 IBM RXN 预测的路线中提取反应 SMARTS，纳入你的自建模板库
4. **记录差距**：注意 IBM RXN 对糖基化反应的处理方式，这是橙皮苷合成的关键
