#!/usr/bin/env python3
"""
运行跨黄酮靶标面板的消融实验。

基于 config/flavonoid_target_panel.csv 中定义的靶标，
对每个靶标运行指定的消融实验。

用法:
    python scripts/run_panel_ablation_experiments.py [--targets hesperidin,naringin] [--experiments 8,9]
"""
import argparse
import csv
import json
import sys
import time
from datetime import datetime
from pathlib import Path

from aizynthfinder.aizynthfinder import AiZynthFinder
from rdkit import Chem
from rdkit.Chem import inchi

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_DIR = PROJECT_ROOT / "config"
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "panel_ablation"
PANEL_FILE = CONFIG_DIR / "flavonoid_target_panel.csv"

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# 实验配置（与 run_ablation_experiments.py 一致）
EXPERIMENTS = {
    8: {
        "name": "flavonoid_virtual_bridge",
        "config": "ablation_flavonoid_virtual_bridge.yml",
        "description": "USPTO + custom flavonoid templates + virtual bridge stock",
        "expansion": ["uspto", "ringbreaker", "flavonoid_reaction_families"],
        "stock": ["virtual_bridge", "trusted_intermediate", "strict_buyable"],
    },
    9: {
        "name": "custom_only_virtual_bridge",
        "config": "ablation_custom_only_virtual_bridge.yml",
        "description": "Custom flavonoid templates only + virtual bridge stock",
        "expansion": ["flavonoid_reaction_families"],
        "stock": ["virtual_bridge", "trusted_intermediate", "strict_buyable"],
    },
}


def load_target_panel() -> list:
    """加载靶标面板配置。"""
    if not PANEL_FILE.exists():
        print(f"错误: 靶标面板文件不存在: {PANEL_FILE}")
        sys.exit(1)

    targets = []
    with open(PANEL_FILE, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            enabled = row.get("enabled", "true").strip().lower()
            if enabled not in ("true", "1", "yes", "y"):
                continue
            if row.get("structure_status") != "verified_pubchem":
                print(f"跳过未验证靶标结构: {row.get('target_name', '')}")
                continue
            mol = Chem.MolFromSmiles(row.get("smiles", ""))
            if mol is None:
                print(f"跳过无效SMILES靶标: {row.get('target_name', '')}")
                continue
            expected_inchikey = row.get("pubchem_inchikey", "")
            actual_inchikey = inchi.MolToInchiKey(mol)
            if expected_inchikey and actual_inchikey != expected_inchikey:
                print(
                    f"跳过InChIKey不匹配靶标: {row.get('target_name', '')} "
                    f"expected={expected_inchikey} actual={actual_inchikey}"
                )
                continue
            targets.append(row)
    return targets


def run_experiment_for_target(exp_id: int, target: dict) -> dict:
    """为单个靶标运行消融实验。"""
    exp = EXPERIMENTS[exp_id]
    config_file = CONFIG_DIR / exp["config"]
    target_name = target['target_name']
    target_smiles = target['smiles']

    print(f"\n{'='*60}")
    print(f"实验 {exp_id}: {exp['name']}")
    print(f"靶标: {target_name}")
    print(f"描述: {exp['description']}")
    print(f"配置: {config_file}")
    print(f"开始: {datetime.now().strftime('%H:%M:%S')}")
    print(f"{'='*60}")

    if not config_file.exists():
        print(f"错误: 配置文件不存在: {config_file}")
        return {"exp_id": exp_id, "target": target_name, "error": "config not found"}

    # 初始化
    finder = AiZynthFinder(configfile=str(config_file))

    # 选择 expansion policy
    available_expansion = list(finder.expansion_policy.items)
    selected_expansion = [e for e in exp["expansion"] if e in available_expansion]
    if not selected_expansion:
        print(f"错误: 无可用的 expansion policy。可用: {available_expansion}")
        return {"exp_id": exp_id, "target": target_name, "error": "no expansion policy"}
    finder.expansion_policy.select(selected_expansion)
    print(f"Expansion: {finder.expansion_policy.selection}")

    # 选择 stock
    available_stock = list(finder.stock.items)
    selected_stock = [s for s in exp["stock"] if s in available_stock]
    if not selected_stock:
        print(f"错误: 无可用的 stock。可用: {available_stock}")
        return {"exp_id": exp_id, "target": target_name, "error": "no stock"}
    finder.stock.select(selected_stock)
    print(f"Stock: {finder.stock.selection}")

    # 设置目标
    finder.target_smiles = target_smiles

    # 运行搜索
    start_time = time.time()
    finder.prepare_tree()
    finder.tree_search()
    finder.build_routes()
    elapsed = time.time() - start_time

    # 提取路线数据
    routes_data = []
    for i, route in enumerate(finder.routes):
        reaction_tree = route['reaction_tree']
        route_dict = reaction_tree.to_dict()
        route_dict['score'] = route.get('score', {})
        route_dict['is_solved'] = route.get('route_metadata', {}).get('is_solved', False)
        route_dict['_policy_name'] = route.get('route_metadata', {}).get('policy_name', '')
        routes_data.append(route_dict)

    # 保存 JSON
    target_dir = OUTPUT_DIR / target_name
    target_dir.mkdir(parents=True, exist_ok=True)
    output_file = target_dir / f"{exp['name']}.json"
    with open(output_file, 'w') as f:
        json.dump(routes_data, f, indent=2)

    # 统计
    n_routes = len(routes_data)
    n_solved = sum(1 for r in routes_data if r.get('is_solved', False))

    print(f"\n完成! 耗时: {elapsed:.1f}s")
    print(f"路线数: {n_routes}, 已解决: {n_solved}")
    print(f"输出: {output_file}")

    return {
        "exp_id": exp_id,
        "target": target_name,
        "target_class": target.get('target_class', ''),
        "expected_sugar_family": target.get('expected_sugar_family', ''),
        "name": exp["name"],
        "description": exp["description"],
        "n_routes": n_routes,
        "n_solved": n_solved,
        "elapsed_s": round(elapsed, 1),
        "output_file": str(output_file),
    }


def main():
    parser = argparse.ArgumentParser(description="运行跨黄酮靶标面板消融实验")
    parser.add_argument(
        "--targets", type=str, default="all",
        help="要运行的靶标名称，逗号分隔（默认全部）"
    )
    parser.add_argument(
        "--experiments", type=str, default="8,9",
        help="要运行的实验编号，逗号分隔（默认 8,9）"
    )
    args = parser.parse_args()

    # 加载靶标面板
    targets = load_target_panel()
    if args.targets != "all":
        target_names = [t.strip() for t in args.targets.split(",")]
        targets = [t for t in targets if t['target_name'] in target_names]

    print(f"靶标数量: {len(targets)}")
    for t in targets:
        print(f"  - {t['target_name']} ({t['target_class']})")

    exp_ids = [int(x.strip()) for x in args.experiments.split(",")]
    print(f"实验编号: {exp_ids}")

    # 运行实验
    all_results = []
    for target in targets:
        for exp_id in exp_ids:
            if exp_id not in EXPERIMENTS:
                print(f"警告: 实验 {exp_id} 不存在，跳过")
                continue
            result = run_experiment_for_target(exp_id, target)
            all_results.append(result)

    # 保存汇总
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_file = OUTPUT_DIR / "panel_ablation_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)

    # 打印汇总表
    print(f"\n{'='*80}")
    print("面板消融实验汇总")
    print(f"{'='*80}")
    print(f"{'靶标':<20} {'实验':<30} {'路线数':>6} {'已解决':>6} {'耗时(s)':>8}")
    print("-" * 75)
    for r in all_results:
        if "error" in r:
            print(f"{r['target']:<20} {r['name']:<30} {'ERROR':>6}")
        else:
            print(f"{r['target']:<20} {r['name']:<30} {r['n_routes']:>6} {r['n_solved']:>6} {r['elapsed_s']:>8.1f}")

    print(f"\n汇总文件: {summary_file}")


if __name__ == "__main__":
    main()
