#!/usr/bin/env python3
"""
批量运行消融实验，评估 stock 和 template 的各自贡献。

9 组实验:
  1. baseline_zinc:              USPTO + ZINC
  2. baseline_strict:            USPTO + strict buyable
  3. flavonoid_zinc:             USPTO + custom + ZINC
  4. flavonoid_strict:           USPTO + custom + strict buyable
  5. custom_only_strict:         custom only + strict buyable
  6. flavonoid_trusted:          USPTO + custom + trusted intermediate
  7. custom_only_trusted:        custom only + trusted intermediate
  8. flavonoid_virtual_bridge:   USPTO + custom + virtual bridge
  9. custom_only_virtual_bridge: custom only + virtual bridge

用法:
    python scripts/run_ablation_experiments.py [--experiments 1,2,3,4,5,6,7,8,9]
"""
import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

from aizynthfinder.aizynthfinder import AiZynthFinder

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_DIR = PROJECT_ROOT / "config"
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "ablation"

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# 橙皮苷 SMILES
HESPERIDIN_SMILES = (
    "COc1ccc(C2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)"
    "[C@H](O)[C@@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O"
)

EXPERIMENTS = {
    1: {
        "name": "baseline_zinc",
        "config": "ablation_baseline_zinc.yml",
        "description": "USPTO + ZINC stock",
        "expansion": ["uspto"],
        "stock": ["zinc"],
    },
    2: {
        "name": "baseline_strict",
        "config": "ablation_baseline_strict.yml",
        "description": "USPTO + strict buyable stock",
        "expansion": ["uspto"],
        "stock": ["strict_buyable"],
    },
    3: {
        "name": "flavonoid_zinc",
        "config": "ablation_flavonoid_zinc.yml",
        "description": "USPTO + mapped flavonoid reaction families + ZINC stock",
        "expansion": ["uspto", "ringbreaker", "flavonoid_reaction_families"],
        "stock": ["zinc"],
    },
    4: {
        "name": "flavonoid_strict",
        "config": "ablation_flavonoid_strict.yml",
        "description": "USPTO + mapped flavonoid reaction families + strict buyable stock",
        "expansion": ["uspto", "ringbreaker", "flavonoid_reaction_families"],
        "stock": ["strict_buyable"],
    },
    5: {
        "name": "custom_only_strict",
        "config": "ablation_custom_only_strict.yml",
        "description": "Mapped flavonoid reaction families only + strict buyable stock",
        "expansion": ["flavonoid_reaction_families"],
        "stock": ["strict_buyable"],
    },
    6: {
        "name": "flavonoid_trusted",
        "config": "ablation_flavonoid_trusted.yml",
        "description": "USPTO + custom flavonoid templates + trusted intermediate stock",
        "expansion": ["uspto", "ringbreaker", "flavonoid_reaction_families"],
        "stock": ["trusted_intermediate", "strict_buyable"],
    },
    7: {
        "name": "custom_only_trusted",
        "config": "ablation_custom_only_trusted.yml",
        "description": "Custom flavonoid templates only + trusted intermediate stock",
        "expansion": ["flavonoid_reaction_families"],
        "stock": ["trusted_intermediate", "strict_buyable"],
    },
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


def run_experiment(exp_id: int) -> dict:
    """运行单个消融实验，返回结果摘要。"""
    exp = EXPERIMENTS[exp_id]
    config_file = CONFIG_DIR / exp["config"]

    print(f"\n{'='*60}")
    print(f"实验 {exp_id}: {exp['name']}")
    print(f"描述: {exp['description']}")
    print(f"配置: {config_file}")
    print(f"开始: {datetime.now().strftime('%H:%M:%S')}")
    print(f"{'='*60}")

    if not config_file.exists():
        print(f"错误: 配置文件不存在: {config_file}")
        return {"exp_id": exp_id, "name": exp["name"], "error": "config not found"}

    # 初始化
    finder = AiZynthFinder(configfile=str(config_file))

    # 选择 expansion policy
    available_expansion = list(finder.expansion_policy.items)
    selected_expansion = [e for e in exp["expansion"] if e in available_expansion]
    if not selected_expansion:
        print(f"错误: 无可用的 expansion policy。可用: {available_expansion}")
        return {"exp_id": exp_id, "name": exp["name"], "error": "no expansion policy"}
    finder.expansion_policy.select(selected_expansion)
    print(f"Expansion: {finder.expansion_policy.selection}")

    # 选择 stock
    available_stock = list(finder.stock.items)
    selected_stock = [s for s in exp["stock"] if s in available_stock]
    if not selected_stock:
        print(f"错误: 无可用的 stock。可用: {available_stock}")
        return {"exp_id": exp_id, "name": exp["name"], "error": "no stock"}
    finder.stock.select(selected_stock)
    print(f"Stock: {finder.stock.selection}")

    # 设置目标
    finder.target_smiles = HESPERIDIN_SMILES

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
        # 保存 policy 信息
        route_dict['_policy_name'] = route.get('route_metadata', {}).get('policy_name', '')
        routes_data.append(route_dict)

    # 保存 JSON
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file = OUTPUT_DIR / f"{exp['name']}.json"
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
        "name": exp["name"],
        "description": exp["description"],
        "n_routes": n_routes,
        "n_solved": n_solved,
        "elapsed_s": round(elapsed, 1),
        "output_file": str(output_file),
    }


def main():
    parser = argparse.ArgumentParser(description="运行消融实验")
    parser.add_argument(
        "--experiments", type=str, default="1,2,3,4,5,6,7,8,9",
        help="要运行的实验编号，逗号分隔（默认全部）"
    )
    args = parser.parse_args()

    exp_ids = [int(x.strip()) for x in args.experiments.split(",")]
    print(f"将运行实验: {exp_ids}")
    print(f"目标分子: {HESPERIDIN_SMILES[:60]}...")

    results = []
    for exp_id in exp_ids:
        if exp_id not in EXPERIMENTS:
            print(f"警告: 实验 {exp_id} 不存在，跳过")
            continue
        result = run_experiment(exp_id)
        results.append(result)

    # Save summary. When running a subset, preserve previous entries for
    # experiments that were not rerun in this invocation.
    summary_file = OUTPUT_DIR / "ablation_summary.json"
    merged_results = {}
    if summary_file.exists():
        try:
            with open(summary_file) as f:
                for item in json.load(f):
                    if "exp_id" in item:
                        merged_results[item["exp_id"]] = item
        except json.JSONDecodeError:
            merged_results = {}
    for item in results:
        if "exp_id" in item:
            merged_results[item["exp_id"]] = item

    summary_results = [merged_results[idx] for idx in sorted(merged_results)]
    with open(summary_file, 'w') as f:
        json.dump(summary_results, f, indent=2, ensure_ascii=False)

    # 打印汇总表
    print(f"\n{'='*60}")
    print("消融实验汇总")
    print(f"{'='*60}")
    print(f"{'实验':<25} {'路线数':>6} {'已解决':>6} {'耗时(s)':>8}")
    print("-" * 50)
    for r in results:
        if "error" in r:
            print(f"{r['name']:<25} {'ERROR':>6}")
        else:
            print(f"{r['name']:<25} {r['n_routes']:>6} {r['n_solved']:>6} {r['elapsed_s']:>8.1f}")

    print(f"\n汇总文件: {summary_file}")


if __name__ == "__main__":
    main()
