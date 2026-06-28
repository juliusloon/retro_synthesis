#!/home/ljx/retro_synthesis/venv/bin/python3
"""
橙皮苷逆合成结果可视化分析脚本（Jupyter Notebook 专用）

使用方法：
1. 确保已安装依赖: pip install rdkit matplotlib seaborn pandas
2. 在 notebook cell 中运行: %run hesperidin_result_visualization.py
3. 或在 notebook 中逐段复制粘贴

功能：
- 路线统计概览（步数、前体分子量、评分分布）
- 单条路线化学结构可视化
- 多条路线对比图表
- 逆合成树层次可视化
- 起始原料商业可得性分析
"""

import json
import os
from collections import Counter

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

# ============================
# 配置
# ============================
RESULT_FILE = "outputs/hesperidin_result.json"
OUTPUT_DIR = "outputs/analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 设置中文字体（如果系统支持）和图表风格
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

# ============================
# 加载数据
# ============================

with open(RESULT_FILE, "r", encoding="utf-8") as f:
    data = json.load(f)

routes = data.get("routes", [])
print(f"✓ 加载 {len(routes)} 条路线")
print(f"  目标分子量: {data['target_mw']}")
print(f"  搜索参数: max_transforms={data['search_parameters']['max_transforms']}, "
      f"iterations={data['search_parameters']['iteration_limit']}")

# ============================
# 1. 路线统计概览（大图表）
# ============================

fig = plt.figure(figsize=(16, 10))
gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

# --- 1a. 步数分布 ---
ax1 = fig.add_subplot(gs[0, 0])
steps = [r["steps"] for r in routes]
step_counts = Counter(steps)
ax1.bar(step_counts.keys(), step_counts.values(), color="steelblue", edgecolor="black")
ax1.set_xlabel("Number of Steps")
ax1.set_ylabel("Number of Routes")
ax1.set_title("Route Length Distribution")
ax1.set_xticks(range(min(steps), max(steps)+1))
for k, v in step_counts.items():
    ax1.text(k, v + 0.05, str(v), ha="center", va="bottom", fontsize=10)

# --- 1b. 最大前体分子量分布 ---
ax2 = fig.add_subplot(gs[0, 1])
precursor_mws = [r["max_precursor_mw"] for r in routes]
colors = ["green" if mw < 300 else "orange" if mw < 500 else "red" for mw in precursor_mws]
ax2.scatter(range(len(precursor_mws)), sorted(precursor_mws, reverse=True),
            c=[colors[i] for i in range(len(precursor_mws))], s=60, edgecolors="black")
ax2.axhline(y=300, color="green", linestyle="--", alpha=0.5, label="Ideal (<300)")
ax2.axhline(y=500, color="orange", linestyle="--", alpha=0.5, label="Acceptable (<500)")
ax2.set_xlabel("Route Rank")
ax2.set_ylabel("Max Precursor MW")
ax2.set_title("Precursor Complexity per Route")
ax2.legend()

# --- 1c. 评分分布 ---
ax3 = fig.add_subplot(gs[0, 2])
scores = [r["score"] for r in routes]
bins = 10 if len(scores) >= 10 else len(scores)
ax3.hist(scores, bins=bins, color="coral", edgecolor="black", alpha=0.7)
ax3.set_xlabel("Route Score")
ax3.set_ylabel("Frequency")
ax3.set_title("Score Distribution")
ax3.axvline(x=max(scores), color="darkgreen", linestyle="--", label=f"Best: {max(scores):.1f}")
ax3.legend()

# --- 1d. 步数 vs 最大前体分子量 散点图 ---
ax4 = fig.add_subplot(gs[1, 0])
scatter = ax4.scatter(steps, precursor_mws, c=scores, cmap="RdYlGn",
                       s=100, edgecolors="black", alpha=0.8)
plt.colorbar(scatter, ax=ax4, label="Route Score")
ax4.set_xlabel("Number of Steps")
ax4.set_ylabel("Max Precursor MW")
ax4.set_title("Steps vs Precursor Complexity")
ax4.axhline(y=300, color="green", linestyle="--", alpha=0.3)
ax4.axhline(y=500, color="orange", linestyle="--", alpha=0.3)

# --- 1e. Stock 可得性饼图 ---
ax5 = fig.add_subplot(gs[1, 1])
stock_counts = Counter([r["all_precursors_in_stock"] for r in routes])
labels = ["All in Stock", "Missing Precursors"]
sizes = [stock_counts.get(True, 0), stock_counts.get(False, 0)]
colors_pie = ["#2ecc71", "#e74c3c"]
ax5.pie(sizes, labels=labels, autopct="%1.0f%%", colors=colors_pie,
        startangle=90, textprops={"fontsize": 11})
ax5.set_title("Stock Availability")

# --- 1f. Top 5 路线评分柱状图 ---
ax6 = fig.add_subplot(gs[1, 2])
top5 = routes[:5]
bars = ax6.barh(range(len(top5)), [r["score"] for r in top5],
                color=["gold", "silver", "#cd7f32", "steelblue", "steelblue"],
                edgecolor="black")
ax6.set_yticks(range(len(top5)))
ax6.set_yticklabels([f"Route {r['route_index']}" for r in top5])
ax6.invert_yaxis()
ax6.set_xlabel("Score")
ax6.set_title("Top 5 Route Scores")
for i, (bar, r) in enumerate(zip(bars, top5)):
    ax6.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
             f"{r['score']:.1f}", va="center", fontsize=10)

fig.suptitle(f"Hesperidin Retrosynthesis Overview ({len(routes)} routes)",
             fontsize=14, fontweight="bold", y=0.98)

plt.savefig(f"{OUTPUT_DIR}/01_overview_statistics.png", dpi=200, bbox_inches="tight")
print(f"✓ 统计概览图已保存: {OUTPUT_DIR}/01_overview_statistics.png")
plt.show()

# ============================
# 2. Top 3 路线详细对比表
# ============================

print("\n" + "="*60)
print("Top 3 路线详细对比")
print("="*60)

fig, axes = plt.subplots(3, 1, figsize=(14, 10))

for idx, route in enumerate(routes[:3]):
    ax = axes[idx]
    ax.axis("off")

    # 构建表格数据
    table_data = []
    headers = ["Step", "Reactants", "→", "Product", "In Stock"]

    for step in route["reactions"]:
        reactants = "\n".join(step["reactants"])
        stock_flags = "\n".join(["✓" if s else "✗" for s in step["in_stock"]])
        table_data.append([
            f"Step {step['step']}",
            reactants[:60] + ("..." if len(reactants) > 60 else ""),
            "→",
            step["product"][:60] + ("..." if len(step["product"]) > 60 else ""),
            stock_flags
        ])

    title = (f"Route {route['route_index']} | Score: {route['score']:.1f} | "
             f"Steps: {route['steps']} | Max Precursor MW: {route['max_precursor_mw']:.0f} | "
             f"All in Stock: {route['all_precursors_in_stock']}")

    table = ax.table(cellText=table_data, colLabels=headers,
                     cellLoc="left", loc="center",
                     colColours=["#4472C4"]*5)
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 2.5)

    # 设置表头样式
    for i in range(5):
        table[(0, i)].set_text_props(color="white", fontweight="bold")
        table[(0, i)].set_facecolor("#4472C4")

    # 交替行颜色
    for i in range(1, len(table_data)+1):
        color = "#E7E6E6" if i % 2 == 0 else "white"
        for j in range(5):
            table[(i, j)].set_facecolor(color)

    ax.set_title(title, fontsize=11, fontweight="bold", pad=10)

plt.suptitle("Top 3 Route Comparison", fontsize=13, fontweight="bold", y=0.98)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/02_top3_route_comparison.png", dpi=200, bbox_inches="tight")
print(f"✓ 路线对比表已保存: {OUTPUT_DIR}/02_top3_route_comparison.png")
plt.show()

# ============================
# 3. Top 1 路线的化学结构可视化
# ============================

print("\n" + "="*60)
print("Top 1 路线化学结构可视化")
print("="*60)

if routes:
    top_route = routes[0]

    fig = plt.figure(figsize=(18, 6))

    all_smiles = []
    all_labels = []

    for step in top_route["reactions"]:
        # 反应物
        for j, smi in enumerate(step["reactants"]):
            all_smiles.append(smi)
            all_labels.append(f"Step {step['step']} Reactant {j+1}")
        # 产物（每步的产物）
        all_smiles.append(step["product"])
        all_labels.append(f"Step {step['step']} Product")

    # 去重但保留顺序
    seen = set()
    unique_smiles = []
    unique_labels = []
    for smi, lab in zip(all_smiles, all_labels):
        if smi not in seen:
            seen.add(smi)
            unique_smiles.append(smi)
            unique_labels.append(lab)

    mols = []
    valid_labels = []
    for smi, lab in zip(unique_smiles, unique_labels):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mols.append(mol)
            valid_labels.append(lab)

    ncols = min(6, len(mols))
    nrows = (len(mols) + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3.5*nrows))
    if nrows == 1 and ncols == 1:
        axes = [[axes]]
    elif nrows == 1 or ncols == 1:
        axes = [axes] if nrows == 1 else [[a] for a in axes]
    else:
        axes = axes.reshape(nrows, ncols)

    for idx, (mol, lab) in enumerate(zip(mols, valid_labels)):
        row = idx // ncols
        col = idx % ncols
        ax = axes[row][col]

        img = Draw.MolToImage(mol, size=(400, 400))
        ax.imshow(img)
        ax.set_title(lab, fontsize=9)
        ax.axis("off")

        # 添加分子量标注
        mw = Descriptors.ExactMolWt(mol)
        ax.text(0.5, -0.05, f"MW: {mw:.1f}", transform=ax.transAxes,
                ha="center", fontsize=8, color="navy")

    # 隐藏多余的子图
    for idx in range(len(mols), nrows * ncols):
        row = idx // ncols
        col = idx % ncols
        axes[row][col].axis("off")

    fig.suptitle(f"Top Route ({top_route['steps']} steps) - Molecular Structures",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/03_top1_molecular_structures.png", dpi=200, bbox_inches="tight")
    print(f"✓ Top 1 分子结构图已保存: {OUTPUT_DIR}/03_top1_molecular_structures.png")
    plt.show()

# ============================
# 4. 逆合成树层次可视化（文本 + 图）
# ============================

print("\n" + "="*60)
print("逆合成树层次结构")
print("="*60)

fig, axes = plt.subplots(min(3, len(routes)), 1, figsize=(14, 4*min(3, len(routes))))
if len(routes) == 1:
    axes = [axes]

for route_idx, route in enumerate(routes[:3]):
    ax = axes[route_idx]
    ax.axis("off")

    # 构建层次文本
    lines = []
    lines.append(f"Route {route['route_index']} (Score: {route['score']:.1f})")
    lines.append("=" * 50)

    for step in reversed(route["reactions"]):  # 从最终产物往回展示
        indent = "  " * (route["steps"] - step["step"])
        reactants_str = " + ".join([f"[{smi[:30]}...]" if len(smi) > 30 else f"[{smi}]"
                                     for smi in step["reactants"]])
        product_str = step["product"][:40] + "..." if len(step["product"]) > 40 else step["product"]
        lines.append(f"{indent}({step['step']}) {reactants_str}")
        lines.append(f"{indent}      → {product_str}")

    text = "\n".join(lines)
    ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=8,
            verticalalignment="top", fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.8))

plt.suptitle("Retrosynthetic Tree Structure (Top 3 Routes)", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/04_tree_structure.png", dpi=200, bbox_inches="tight")
print(f"✓ 逆合成树结构图已保存: {OUTPUT_DIR}/04_tree_structure.png")
plt.show()

# ============================
# 5. 起始原料分析
# ============================

print("\n" + "="*60)
print("起始原料分析")
print("="*60)

all_precursors = []
for route in routes:
    for step in route["reactions"]:
        for j, smi in enumerate(step["reactants"]):
            if step["in_stock"][j]:
                mol = Chem.MolFromSmiles(smi)
                mw = Descriptors.ExactMolWt(mol) if mol else 0
                all_precursors.append({
                    "smiles": smi,
                    "mw": mw,
                    "route_idx": route["route_index"],
                    "step": step["step"]
                })

# 去重统计
unique_precursors = {}
for p in all_precursors:
    smi = p["smiles"]
    if smi not in unique_precursors:
        unique_precursors[smi] = {"mw": p["mw"], "count": 0, "routes": set()}
    unique_precursors[smi]["count"] += 1
    unique_precursors[smi]["routes"].add(p["route_idx"])

print(f"\n共发现 {len(unique_precursors)} 种独特起始原料")
print("\n按出现频率排序的前 10 种:")

sorted_precursors = sorted(unique_precursors.items(), key=lambda x: x[1]["count"], reverse=True)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# 左图：前 10 种原料的出现频率
top10 = sorted_precursors[:10]
names = [f"MW: {p[1]['mw']:.0f}" for p in top10]
counts = [p[1]["count"] for p in top10]
colors = plt.cm.viridis([0.1 + 0.8 * i / len(top10) for i in range(len(top10))])

bars = ax1.barh(range(len(top10)), counts, color=colors, edgecolor="black")
ax1.set_yticks(range(len(top10)))
ax1.set_yticklabels(names)
ax1.invert_yaxis()
ax1.set_xlabel("Frequency (appears in N routes)")
ax1.set_title("Top 10 Starting Materials by Frequency")

for i, (bar, (smi, info)) in enumerate(zip(bars, top10)):
    # 在柱子右侧标注 SMILES 缩写
    short_smi = smi[:40] + "..." if len(smi) > 40 else smi
    ax1.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
             short_smi, va="center", fontsize=7)

# 右图：分子量分布
all_mws = [p["mw"] for p in unique_precursors.values()]
ax2.hist(all_mws, bins=15, color="mediumpurple", edgecolor="black", alpha=0.7)
ax2.set_xlabel("Molecular Weight")
ax2.set_ylabel("Number of Precursors")
ax2.set_title("Starting Material MW Distribution")
ax2.axvline(x=300, color="green", linestyle="--", label="Ideal (<300)")
ax2.axvline(x=sum(all_mws)/len(all_mws), color="red", linestyle="--",
            label=f"Mean: {sum(all_mws)/len(all_mws):.0f}")
ax2.legend()

plt.suptitle("Starting Material Analysis", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/05_starting_materials.png", dpi=200, bbox_inches="tight")
print(f"✓ 起始原料分析图已保存: {OUTPUT_DIR}/05_starting_materials.png")
plt.show()

# ============================
# 6. 生成 Markdown 分析报告
# ============================

print("\n" + "="*60)
print("生成 Markdown 分析报告")
print("="*60)

md_content = f"""# 橙皮苷逆合成分析报告

## 搜索参数

| 参数 | 值 |
|------|-----|
| 目标分子 | Hesperidin |
| 分子量 | {data['target_mw']:.1f} |
| 搜索深度 (max_transforms) | {data['search_parameters']['max_transforms']} |
| 迭代次数 | {data['search_parameters']['iteration_limit']} |
| 时间限制 | {data['search_parameters']['time_limit']}s |
| UCB 常数 C | {data['search_parameters']['C']} |
| 模板数 cutoff | {data['search_parameters']['cutoff_number']} |

## 搜索结果统计

| 指标 | 值 |
|------|-----|
| 搜索耗时 | {data['search_stats']['elapsed_time_sec']:.1f}s |
| 原始路线数 | {data['search_stats']['total_routes_found']} |
| 过滤后路线数 | {data['search_stats']['filtered_routes']} |
| 最佳路线评分 | {routes[0]['score'] if routes else 'N/A':.1f} |

## Top 5 路线摘要

"""

for i, route in enumerate(routes[:5], 1):
    md_content += f"""### Route {route['route_index']} (Rank {i})

- **评分**: {route['score']:.1f}
- **步数**: {route['steps']}
- **最大前体分子量**: {route['max_precursor_mw']:.1f}
- **所有前体在 stock 中**: {'✓' if route['all_precursors_in_stock'] else '✗'}

| Step | Reactants | Product | In Stock |
|------|-----------|---------|----------|
"""
    for step in route["reactions"]:
        reactants = "<br>".join(step["reactants"])
        stock = "<br>".join(["✓" if s else "✗" for s in step["in_stock"]])
        md_content += f"| {step['step']} | {reactants} | {step['product']} | {stock} |\n"
    md_content += "\n"

md_content += """## 可视化图表

- [统计概览](01_overview_statistics.png)
- [路线对比](02_top3_route_comparison.png)
- [分子结构](03_top1_molecular_structures.png)
- [树结构](04_tree_structure.png)
- [起始原料](05_starting_materials.png)

---

*Generated by hesperidin_result_visualization.py*
"""

md_path = f"{OUTPUT_DIR}/analysis_report.md"
with open(md_path, "w", encoding="utf-8") as f:
    f.write(md_content)

print(f"✓ Markdown 报告已保存: {md_path}")

print("\n" + "="*60)
print("✓ 所有可视化分析完成！")
print(f"  输出目录: {OUTPUT_DIR}/")
print("  生成的文件:")
for fname in sorted(os.listdir(OUTPUT_DIR)):
    print(f"    - {fname}")
print("="*60)