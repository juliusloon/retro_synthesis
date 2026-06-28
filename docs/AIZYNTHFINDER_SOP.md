# AiZynthFinder 终端运行标准操作程序（SOP）

> **适用环境**：Ubuntu Linux（本地 / WSL / 服务器）  
> **项目目录**：`/home/ljx/retro_synthesis`  
> **Python 版本**：3.11（推荐，官方支持 3.10–3.12）  
> **更新日期**：2026-06-04

---

## 1. 环境准备

### 1.1 激活虚拟环境

所有操作均在项目目录下执行。先进入目录并激活已创建的虚拟环境：

```bash
cd /home/ljx/retro_synthesis
source venv/bin/activate
```

激活后，命令行提示符前会出现 `(venv)` 标识：

```
(venv) ljx@hostname:~/retro_synthesis$
```

### 1.2 安装 AiZynthFinder

> ⚠️ **Python 版本兼容性警告**  
> 当前 `venv` 默认 Python 为 3.13，但 **AiZynthFinder 官方仅支持 3.10–3.12**。直接 `pip install` 会报 `No matching distribution` 错误。  
> **解决方案**：使用项目已安装的 Miniforge3 新建一个 Python 3.11 环境，并在该环境中安装：
>
> ```bash
> /home/ljx/miniforge3/bin/conda create -n retro python=3.11 -y
> /home/ljx/miniforge3/bin/conda run -n retro pip install aizynthfinder[all]
> # 后续激活方式改为：
> /home/ljx/miniforge3/bin/conda activate retro
> ```
> 下文示例仍使用原有 `venv`，若您已按上述步骤创建 `retro` conda 环境，请将所有 `source venv/bin/activate` 替换为 `/home/ljx/miniforge3/bin/conda activate retro`。

在正确版本的 Python 环境中执行：

```bash
# 确保在 venv / conda 激活状态下执行
pip install --upgrade pip
pip install aizynthfinder[all]
```

安装完成后验证：

```bash
aizynthcli --version
```

若显示版本号，则安装成功。

---

## 2. 项目文件结构说明

本项目已预置模型与配置文件，目录结构如下：

```
retro_synthesis/
├── venv/                          # Python 虚拟环境
├── config/
│   ├── data.yml                   # 主配置文件（已配置好模型路径）
│   ├── flavonoid.yml              # 黄酮专用配置
│   └── baseline.yml               # 基线配置
├── data/
│   ├── uspto_model.onnx           # USPTO 扩展策略模型
│   ├── uspto_templates.csv.gz     # USPTO 反应模板
│   ├── uspto_ringbreaker_model.onnx      # Ringbreaker 策略模型
│   ├── uspto_ringbreaker_templates.csv.gz # Ringbreaker 模板
│   ├── uspto_filter_model.onnx    # 过滤器模型
│   ├── zinc_stock_17_04_20.hdf5   # ZINC 库存数据
│   └── zinc_stock.hdf5 -> ...     # 库存软链接
├── src/                           # 源代码（模板提取/验证流水线）
├── scripts/                       # 独立脚本（custom_expansion, processing）
├── templates/                     # 自定义黄酮反应模板库
├── notebooks/                     # Jupyter notebooks
├── outputs/                       # 搜索结果
└── docs/                          # 文档
```

### 2.1 配置文件 `config/data.yml`

当前已配置的内容如下，**通常无需修改**：

```yaml
expansion:
  uspto:
    - /home/ljx/retro_synthesis/data/uspto_model.onnx
    - /home/ljx/retro_synthesis/data/uspto_templates.csv.gz
  ringbreaker:
    - /home/ljx/retro_synthesis/data/uspto_ringbreaker_model.onnx
    - /home/ljx/retro_synthesis/data/uspto_ringbreaker_templates.csv.gz
filter:
  uspto: /home/ljx/retro_synthesis/data/uspto_filter_model.onnx
stock:
  zinc: /home/ljx/retro_synthesis/data/zinc_stock.hdf5
```

| 字段 | 说明 |
|------|------|
| `expansion` | 定义扩展策略（神经网络模型 + 反应模板）。`uspto` 为标准策略，`ringbreaker` 用于开环反应。 |
| `filter` | 定义过滤模型，用于剔除低质量反应建议。 |
| `stock` | 定义可购买 / 可获取的起始原料库存。本项目使用 ZINC 数据库子集。 |

---

## 3. 运行方式

### 3.1 方式一：命令行直接输入 SMILES（快速测试）

适合对单个分子进行快速逆合成分析。

```bash
cd /home/ljx/retro_synthesis
source venv/bin/activate

aizynthcli --smiles "COc1ccc([C@@H]2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)[C@H](O)[C@@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O" \
           --config config/data.yml \
           --output output_single.json
```

> 上述 SMILES 来自 `smiles.txt`，对应分子为 **橙皮苷（Hesperidin）**。

### 3.2 方式二：从文件批量输入（推荐）

若需分析多个分子，将 SMILES 写入文本文件（每行一个），然后批量读取。

**示例**：`smiles.txt` 已存在于项目根目录，当前内容为：

```
COc1ccc([C@@H]2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)[C@H](O)[C@@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O
```

运行命令：

```bash
cd /home/ljx/retro_synthesis
source venv/bin/activate

aizynthcli --smiles smiles.txt \
           --config config/data.yml \
           --output output_routes.json
```

若 SMILES 文件有多行，程序会自动逐条解析。

### 3.3 方式三：使用 SDF/MOL 文件输入

对于已保存为结构文件的分子（如 `Hesperidin.sdf`）：

```bash
aizynthcli --smiles Hesperidin.sdf \
           --config config/data.yml \
           --output output_sdf.json
```

---

## 4. 常用参数速查

| 参数 | 说明 | 示例 |
|------|------|------|
| `--smiles` | 输入 SMILES 字符串或文件路径 | `--smiles "CCO"` 或 `--smiles input.txt` |
| `--config` | 指定配置文件路径 | `--config config/data.yml` |
| `--output` | 输出 JSON 文件路径 | `--output result.json` |
| `--nclusters` | 返回的逆合成路线簇数量 | `--nclusters 5`（默认） |
| `--C` | 蒙特卡洛 搜索次数（树探索次数） | `--C 100`（默认 100） |
| `--time_limit` | 单分子的搜索时间限制（秒） | `--time_limit 120` |
| `--no_stats` | 不输出统计信息到终端 | `--no_stats` |

**完整示例**（提高搜索强度）：

```bash
aizynthcli --smiles smiles.txt \
           --config config/data.yml \
           --output output_intense.json \
           --C 200 \
           --time_limit 300 \
           --nclusters 10
```

---

## 5. 输出结果解读

运行成功后，会生成一个 JSON 文件（如 `output_routes.json`）。其顶层结构如下：

```json
{
  "request": {
    "smiles": "...",
    "config": "...",
    "created_at": "2026-06-04T..."
  },
  "routes": [
    {
      "cluster_id": 0,
      "score": 0.95,
      "reactions": [...]
    }
  ],
  "stats": {
    "search_time": 45.2,
    "first_solution_time": 12.3
  }
}
```

| 字段 | 含义 |
|------|------|
| `routes` | 找到的逆合成路线列表，按 `score` 排序。 |
| `score` | 路线得分（0~1），越高表示模型认为该路线越可靠。 |
| `cluster_id` | 路线聚类 ID，相同 ID 表示化学反应策略相似。 |
| `reactions` | 该路线包含的具体反应步骤（SMILES 表示）。 |
| `stats.search_time` | 实际搜索耗时（秒）。 |

**可视化路线**：可使用 AiZynthFinder 的 Web 界面（见第 6 章），或提取 JSON 中的 `reactions` 字段用 RDKit 绘制。

---

## 6. Web/GUI 交互式界面（Jupyter Notebook）

AiZynthFinder 官方提供的图形界面 **`AiZynthApp`** 是基于 **Jupyter Notebook** 的交互式应用。相比 `aizynthcli` 的批量输出，GUI 更适合对**单个分子**进行深入的、可视化的逆合成路线筛选与对比。

### 6.1 前置依赖

GUI 需要 `jupyter` 支持，请确保已安装：

```bash
pip install notebook
```

### 6.2 方式一：一键自动创建 Notebook（推荐）

`aizynthfinder` 提供了命令行工具 `aizynthapp`，可自动生成包含 GUI 代码的 Notebook 并打开浏览器：

```bash
cd /home/ljx/retro_synthesis
source venv/bin/activate

aizynthapp --config config/data.yml
```

执行后：
1. 当前目录会自动生成一个 `.ipynb` 文件（如 `aizynth_app.ipynb`）。
2. 系统会尝试自动打开浏览器并跳转到该 Notebook。
3. 在第一个 Cell 中点击 **Run**（或按 `Shift+Enter`），即可看到交互式界面。

> **服务器用户注意**：若程序提示 `No web browser found`，请改用**方式二**手动启动 Jupyter，并通过端口转发访问。

### 6.3 方式二：手动创建并启动 Notebook

若 `aizynthapp` 不可用，或您需要自定义 Jupyter 端口，请按以下步骤操作：

**步骤 1：启动 Jupyter Notebook**

```bash
cd /home/ljx/retro_synthesis
source venv/bin/activate

jupyter notebook --no-browser --ip=0.0.0.0 --port=8888 --NotebookApp.token='retro2024'
```

参数说明：
- `--no-browser`：不在服务器上自动弹窗（无图形界面服务器必需）。
- `--ip=0.0.0.0`：允许外部 IP 访问（配合 SSH 端口转发使用）。
- `--port=8888`：指定服务端口。
- `--NotebookApp.token='retro2024'`：设置访问密码（可选，建议生产环境使用随机 Token）。

**步骤 2：本地端口转发（远程服务器必需）**

在**本地电脑**的终端（Windows PowerShell / macOS Terminal / WSL）执行：

```bash
ssh -L 8888:localhost:8888 ljx@<your_server_ip>
```

然后在本地浏览器访问：

```
http://localhost:8888/?token=retro2024
```

**步骤 3：新建 Notebook 并加载 GUI**

在 Jupyter 页面点击 **New → Python 3 (ipykernel)**，在第一个 Cell 中输入：

```python
from aizynthfinder.interfaces import AiZynthApp
app = AiZynthApp("config/data.yml")
```

按 `Shift+Enter` 运行，稍等片刻即可看到交互式面板。

### 6.4 界面操作流程

GUI 加载后，按以下步骤操作：

1. **输入目标分子**：在 **Target SMILES** 文本框中粘贴 SMILES（例如 `smiles.txt` 中的橙皮苷）。
2. **选择库存（Stocks）**：下拉菜单勾选 `zinc`。
3. **选择策略（Policy）**：下拉菜单勾选 `uspto`；如需处理环状结构，可同时勾选 `ringbreaker`。
4. **调整搜索参数**（可选）：
   - **C**：蒙特卡洛探索次数（默认 100，可调至 200 提高覆盖率）。
   - **Time limit**：单分子搜索时间上限（秒）。
   - **Max iterations**：最大迭代次数。
5. 点击 **Run Search** 开始树搜索。
6. 搜索完成后：
   - 点击 **Show Reactions** 查看排名最高的逆合成路线列表。
   - 点击 **Show Tree** 查看完整的蒙特卡洛搜索树可视化。
   - 使用 **Sort by** 下拉框可切换不同的评分函数对路线重新排序。

### 6.5 结果深入分析（Python API）

在 GUI 后续的 Cell 中，可通过 `app.finder` 访问底层 `AiZynthFinder` 对象，提取更详细的统计信息：

```python
finder = app.finder
stats = finder.extract_statistics()
print(stats)
```

获取所有找到的路线：

```python
routes = finder.routes
print(f"共找到 {len(routes)} 条路线")
for i, route in enumerate(routes[:5]):
    print(f"Route {i+1}: score = {route.score}")
```

### 6.6 路线聚类分析

若路线数量较多，可使用内置的聚类工具分析路线间的相似性：

在新 Cell 中输入：

```python
%matplotlib inline
from aizynthfinder.interfaces.gui.clustering import ClusteringGui
ClusteringGui.from_app(app)
```

运行后会显示一个层次聚类（Dendrogram）交互界面：
- 拖动 **N clusters** 滑块选择聚类数量。
- 系统会将策略相似的路线归为一组，便于您对比不同化学策略的优劣。

### 6.7 CLI 与 GUI 对比

| 特性 | `aizynthcli`（终端） | `AiZynthApp`（GUI） |
|------|----------------------|---------------------|
| 适用场景 | 批量分子、自动化脚本 | 单分子深度分析、教学演示 |
| 可视化 | 无（仅 JSON） | 搜索树、反应路线、聚类图 |
| 交互性 | 无 | 可实时调整参数并重新搜索 |
| 输出 | JSON 文件 | Notebook Cell + JSON 导出 |

---

## 7. 退出与清理

### 7.1 退出虚拟环境

```bash
deactivate
```

### 7.2 压缩/归档结果

```bash
cd /home/ljx/retro_synthesis
gzip -k output_routes.json   # 保留原文件并生成 .gz
```

---

## 8. 常见问题排查（FAQ）

### Q1: 提示 `aizynthcli: command not found`
- **原因**：未激活虚拟环境，或 `aizynthfinder` 未安装。
- **解决**：执行 `source venv/bin/activate`，然后 `pip install aizynthfinder[all]`。若 Python 为 3.13，请先创建 Python 3.11 的 conda 环境（见 1.2 节）。

### Q6: `pip install aizynthfinder[all]` 报错 `No matching distribution found`
- **原因**：Python 3.13 不在官方支持范围内。
- **解决**：使用 miniforge3 创建 Python 3.11 环境并重新安装：
  ```bash
  /home/ljx/miniforge3/bin/conda create -n retro python=3.11 -y
  /home/ljx/miniforge3/bin/conda run -n retro pip install aizynthfinder[all]
  ```

### Q7: Jupyter Notebook 中运行 `AiZynthApp` 后无界面显示
- **原因**：部分 Jupyter 前端（如 VS Code 内置 Notebook）对 ipywidgets 支持不完整。
- **解决**：使用传统浏览器打开 Jupyter（`jupyter notebook`），或在 Cell 开头加 `%matplotlib inline` 并确保已安装 `ipywidgets`：
  ```bash
  pip install ipywidgets
  ```

### Q8: 远程服务器上浏览器无法自动打开
- **原因**：服务器无图形界面，`aizynthapp` 无法调用本地浏览器。
- **解决**：使用方式二手动启动 Jupyter，并通过 SSH 端口转发（`-L 8888:localhost:8888`）在本地浏览器访问（见 6.3 节）。

### Q2: 运行时报 `ONNXRuntimeError` 或模型加载失败
- **原因**：模型文件路径配置错误，或文件损坏。
- **解决**：检查 `config/data.yml` 中的绝对路径是否与 `data/` 目录下的实际文件名一致。本项目已使用绝对路径，通常不会出现此问题。

### Q3: 找不到起始原料（`No stock compound`）
- **原因**：目标分子的片段不在 `zinc_stock.hdf5` 库存中。
- **解决**：可尝试更换或扩展库存文件（需自行准备更大的 HDF5 库存），或适当放宽搜索参数。

### Q4: 运行速度极慢 / 内存不足
- **原因**：默认搜索参数对大型分子或低配服务器压力较大。
- **解决**：降低 `--C`（如 50）和 `--time_limit`（如 60），或改用更小规模的库存。

### Q5: `smiles.txt` 中有空行或注释导致解析失败
- **解决**：确保输入文件每行一个有效 SMILES，不要包含空行或 `#` 注释。可用以下命令清理：
  ```bash
  grep -v '^$' smiles.txt | grep -v '^#' > smiles_clean.txt
  ```

---

## 9. 快速开始 checklist

- [ ] `cd /home/ljx/retro_synthesis`
- [ ] `source venv/bin/activate`（或激活 conda `retro` 环境）
- [ ] （首次）`pip install aizynthfinder[all]`（若 Python 3.13 需先建 conda 环境，见 1.2 节）
- [ ] 检查 `config/data.yml` 路径是否正确
- [ ] 准备 SMILES 到 `smiles.txt`
- [ ] （批量分析）执行 `aizynthcli --smiles smiles.txt --config config/data.yml --output output.json`
- [ ] （交互分析）执行 `aizynthapp --config config/data.yml` 或手动启动 Jupyter Notebook
- [ ] 检查生成的 `output.json` 并提取 `routes`，或在 GUI 中查看搜索树

---

## 附录 A：单次完整运行脚本（CLI）

为方便复用，可将以下命令保存为 `run.sh`：

```bash
#!/bin/bash
set -e

cd /home/ljx/retro_synthesis
source venv/bin/activate

CONFIG="config/data.yml"
INPUT="smiles.txt"
OUTPUT="output_$(date +%Y%m%d_%H%M%S).json"

echo "[INFO] 开始逆合成分析..."
echo "[INFO] 输入文件: $INPUT"
echo "[INFO] 配置文件: $CONFIG"
echo "[INFO] 输出文件: $OUTPUT"

aizynthcli --smiles "$INPUT" \
           --config "$CONFIG" \
           --output "$OUTPUT" \
           --C 100 \
           --time_limit 120

echo "[INFO] 分析完成，结果保存至: $OUTPUT"
```

赋予执行权限并运行：

```bash
chmod +x run.sh
./run.sh
```

---

## 附录 B：Jupyter Notebook 快速模板

若选择手动方式启动 GUI，可将以下内容保存为 `aizynth_gui.ipynb` 或直接在 Jupyter Cell 中按顺序执行：

```python
# Cell 1: 加载 GUI
from aizynthfinder.interfaces import AiZynthApp
app = AiZynthApp("config/data.yml")
```

```python
# Cell 2: 搜索完成后提取统计信息（可选）
finder = app.finder
stats = finder.extract_statistics()
print(stats)
```

```python
# Cell 3: 路线聚类分析（可选）
%matplotlib inline
from aizynthfinder.interfaces.gui.clustering import ClusteringGui
ClusteringGui.from_app(app)
```
