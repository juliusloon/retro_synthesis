#!/usr/bin/env python3
"""
比较跨黄酮靶标面板的消融实验结果。

基于 outputs/panel_ablation/ 中的结果文件，
生成靶标感知的比较报告。

用法:
    python scripts/compare_panel_ablation_results.py
"""
import json
import sys
import csv
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "panel_ablation"
STOCK_LAYERS_DIR = PROJECT_ROOT / "templates" / "stock_layers"


def load_stock_layer_inchikeys() -> dict:
    """加载各 stock 层级的 InChIKey 集合。"""
    stock_layers = {}
    for layer_name in ["strict_buyable", "trusted_intermediate", "virtual_bridge"]:
        inchikey_file = STOCK_LAYERS_DIR / f"{layer_name}_stock_inchikeys.txt"
        if inchikey_file.exists():
            with open(inchikey_file) as f:
                stock_layers[layer_name] = set(line.strip() for line in f if line.strip())
        else:
            stock_layers[layer_name] = set()
    return stock_layers


def route_all_leaves_in_layers(leaf_inchikeys: list, stock_layers: dict, layer_names: list) -> bool:
    """检查路线所有叶子是否都属于指定 stock 层级集合。"""
    if not leaf_inchikeys:
        return False
    allowed = set()
    for layer_name in layer_names:
        allowed.update(stock_layers.get(layer_name, set()))
    return all(inchikey in allowed for inchikey in leaf_inchikeys)


def analyze_route(route_dict: dict, stock_layers: dict, exp_name: str) -> dict:
    """分析单条路线。"""
    def collect_leaf_inchikeys(node):
        leaves = []
        if isinstance(node, dict):
            if node.get('type') == 'mol':
                children = node.get('children', [])
                if not children or len(children) == 0:
                    inchikey = node.get('inchikey', '')
                    if inchikey:
                        leaves.append(inchikey)
                else:
                    for child in children:
                        leaves.extend(collect_leaf_inchikeys(child))
            else:
                for child in node.get('children', []):
                    leaves.extend(collect_leaf_inchikeys(child))
        return leaves

    leaf_inchikeys = collect_leaf_inchikeys(route_dict)
    is_solved = route_dict.get('is_solved', False)

    # 检查是否使用虚拟桥接
    uses_virtual_bridge = any(
        ik in stock_layers.get('virtual_bridge', set())
        for ik in leaf_inchikeys
    )

    # 确定路由有效性分类
    if is_solved:
        if uses_virtual_bridge:
            route_validity_class = 'bridge_closed_connectivity'
        elif "zinc" in exp_name:
            route_validity_class = 'zinc_baseline_solved'
        elif route_all_leaves_in_layers(
            leaf_inchikeys,
            stock_layers,
            ["strict_buyable", "trusted_intermediate"],
        ):
            route_validity_class = 'strict_trusted_solved'
        else:
            route_validity_class = 'zinc_baseline_solved'
    else:
        route_validity_class = 'unsolved'

    return {
        'is_solved': is_solved,
        'uses_virtual_bridge': uses_virtual_bridge,
        'route_validity_class': route_validity_class,
        'leaf_inchikeys': leaf_inchikeys,
    }


def load_target_panel() -> dict:
    """加载靶标面板配置，返回靶标名称到配置的映射。"""
    panel_file = PROJECT_ROOT / "config" / "flavonoid_target_panel.csv"
    if not panel_file.exists():
        return {}

    panel = {}
    with open(panel_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            enabled = row.get("enabled", "true").strip().lower()
            if enabled not in ("true", "1", "yes", "y"):
                continue
            if row.get("structure_status") and row.get("structure_status") != "verified_pubchem":
                continue
            panel[row['target_name']] = row
    return panel


def main():
    if not OUTPUT_DIR.exists():
        print("错误: 面板消融结果目录不存在")
        print("请先运行: python scripts/run_panel_ablation_experiments.py")
        sys.exit(1)

    stock_layers = load_stock_layer_inchikeys()
    target_panel = load_target_panel()

    # 收集所有结果
    all_results = {}
    for target_dir in OUTPUT_DIR.iterdir():
        if not target_dir.is_dir():
            continue
        target_name = target_dir.name
        all_results[target_name] = {}

        for result_file in target_dir.glob("*.json"):
            exp_name = result_file.stem
            with open(result_file) as f:
                routes = json.load(f)

            analyses = [analyze_route(r, stock_layers, exp_name) for r in routes]

            n_total = len(routes)
            n_solved = sum(1 for a in analyses if a['is_solved'])
            n_bridge_closed = sum(1 for a in analyses if a['is_solved'] and a['uses_virtual_bridge'])
            n_non_virtual = sum(1 for a in analyses if a['is_solved'] and not a['uses_virtual_bridge'])

            # 路由有效性分类统计
            validity_counts = defaultdict(int)
            for a in analyses:
                validity_counts[a['route_validity_class']] += 1

            all_results[target_name][exp_name] = {
                'n_total': n_total,
                'n_solved': n_solved,
                'n_bridge_closed': n_bridge_closed,
                'n_non_virtual': n_non_virtual,
                'validity_counts': dict(validity_counts),
            }

    # 打印汇总
    print("=" * 100)
    print("跨黄酮靶标面板消融实验结果比较")
    print("=" * 100)

    # 按靶标汇总
    print("\n### 按靶标汇总\n")
    print(f"{'靶标':<20} {'类别':<25} {'糖家族':<15} {'实验':<30} {'解决':>6} {'桥接':>6} {'非虚拟':>6}")
    print("-" * 110)

    for target_name, exps in sorted(all_results.items()):
        target_info = target_panel.get(target_name, {})
        target_class = target_info.get('target_class', 'unknown')
        sugar_family = target_info.get('expected_sugar_family', 'unknown')

        for exp_name, stats in sorted(exps.items()):
            print(f"{target_name:<20} {target_class:<25} {sugar_family:<15} "
                  f"{exp_name:<30} {stats['n_solved']:>6} {stats['n_bridge_closed']:>6} "
                  f"{stats['n_non_virtual']:>6}")

    # 聚合统计
    print("\n### 聚合统计\n")

    # 按实验聚合
    exp_aggregates = defaultdict(lambda: {'n_total': 0, 'n_solved': 0, 'n_bridge_closed': 0, 'n_non_virtual': 0})
    for target_name, exps in all_results.items():
        for exp_name, stats in exps.items():
            exp_aggregates[exp_name]['n_total'] += stats['n_total']
            exp_aggregates[exp_name]['n_solved'] += stats['n_solved']
            exp_aggregates[exp_name]['n_bridge_closed'] += stats['n_bridge_closed']
            exp_aggregates[exp_name]['n_non_virtual'] += stats['n_non_virtual']

    print(f"{'实验':<35} {'总路线':>8} {'解决':>8} {'桥接':>8} {'非虚拟':>8}")
    print("-" * 70)
    for exp_name, agg in sorted(exp_aggregates.items()):
        print(f"{exp_name:<35} {agg['n_total']:>8} {agg['n_solved']:>8} "
              f"{agg['n_bridge_closed']:>8} {agg['n_non_virtual']:>8}")

    # 按糖家族统计
    print("\n### 按糖家族统计\n")
    sugar_family_stats = defaultdict(lambda: {'n_total': 0, 'n_solved': 0, 'n_bridge_closed': 0})
    for target_name, exps in all_results.items():
        target_info = target_panel.get(target_name, {})
        sugar_family = target_info.get('expected_sugar_family', 'unknown')

        for exp_name, stats in exps.items():
            sugar_family_stats[sugar_family]['n_total'] += stats['n_total']
            sugar_family_stats[sugar_family]['n_solved'] += stats['n_solved']
            sugar_family_stats[sugar_family]['n_bridge_closed'] += stats['n_bridge_closed']

    print(f"{'糖家族':<20} {'总路线':>8} {'解决':>8} {'桥接':>8}")
    print("-" * 50)
    for family, stats in sorted(sugar_family_stats.items()):
        print(f"{family:<20} {stats['n_total']:>8} {stats['n_solved']:>8} {stats['n_bridge_closed']:>8}")

    # 保存报告
    report_file = OUTPUT_DIR / "panel_ablation_report.md"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("# 跨黄酮靶标面板消融实验报告\n\n")
        f.write(f"生成日期: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")

        f.write("## 靶标列表\n\n")
        f.write("| 靶标 | 类别 | 预期糖家族 |\n")
        f.write("|------|------|------------|\n")
        for target_name, info in sorted(target_panel.items()):
            f.write(f"| {target_name} | {info.get('target_class', '')} | {info.get('expected_sugar_family', '')} |\n")

        f.write("\n## 详细结果\n\n")
        for target_name, exps in sorted(all_results.items()):
            f.write(f"### {target_name}\n\n")
            f.write(f"- 类别: {target_panel.get(target_name, {}).get('target_class', 'unknown')}\n")
            f.write(f"- 预期糖家族: {target_panel.get(target_name, {}).get('expected_sugar_family', 'unknown')}\n\n")

            f.write("| 实验 | 总路线 | 解决 | 桥接 | 非虚拟 |\n")
            f.write("|------|-------:|-----:|-----:|-------:|\n")
            for exp_name, stats in sorted(exps.items()):
                f.write(f"| {exp_name} | {stats['n_total']} | {stats['n_solved']} | "
                        f"{stats['n_bridge_closed']} | {stats['n_non_virtual']} |\n")
            f.write("\n")

        f.write("## 聚合统计\n\n")
        f.write("### 按实验\n\n")
        f.write("| 实验 | 总路线 | 解决 | 桥接 | 非虚拟 |\n")
        f.write("|------|-------:|-----:|-----:|-------:|\n")
        for exp_name, agg in sorted(exp_aggregates.items()):
            f.write(f"| {exp_name} | {agg['n_total']} | {agg['n_solved']} | "
                    f"{agg['n_bridge_closed']} | {agg['n_non_virtual']} |\n")

        f.write("\n### 按糖家族\n\n")
        f.write("| 糖家族 | 总路线 | 解决 | 桥接 |\n")
        f.write("|--------|-------:|-----:|-----:|\n")
        for family, stats in sorted(sugar_family_stats.items()):
            f.write(f"| {family} | {stats['n_total']} | {stats['n_solved']} | {stats['n_bridge_closed']} |\n")

    # 保存 JSON 报告
    json_report = {
        'generated': __import__('datetime').datetime.now().isoformat(),
        'target_panel': target_panel,
        'results': all_results,
        'exp_aggregates': dict(exp_aggregates),
        'sugar_family_stats': dict(sugar_family_stats),
    }
    json_file = OUTPUT_DIR / "panel_ablation_report.json"
    with open(json_file, 'w') as f:
        json.dump(json_report, f, indent=2, ensure_ascii=False)

    print(f"\n报告已保存: {report_file}")
    print(f"JSON 报告: {json_file}")


if __name__ == "__main__":
    main()
