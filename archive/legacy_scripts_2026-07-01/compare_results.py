#!/usr/bin/env python3
"""
比较baseline和flavonoid配置的逆合成结果
"""
import json
import sys
from pathlib import Path
from collections import defaultdict


def load_routes(filepath: str) -> list:
    """加载路线JSON文件"""
    with open(filepath) as f:
        routes = json.load(f)
    return routes


def analyze_route(route_dict: dict) -> dict:
    """分析单条路线"""
    # 计算步骤数（通过计算reaction节点）
    def count_reactions(node):
        count = 0
        if isinstance(node, dict):
            if 'children' in node:
                count += 1
                for child in node['children']:
                    count += count_reactions(child)
        return count

    # 检查是否完全解决
    def is_solved(node):
        if isinstance(node, dict):
            if 'children' in node:
                if len(node['children']) == 0:
                    return True
                return all(is_solved(child) for child in node['children'])
        return True

    n_steps = count_reactions(route_dict)
    solved = is_solved(route_dict)

    return {
        'n_steps': n_steps,
        'is_solved': solved,
    }


def compare_results(baseline_file: str, flavonoid_file: str):
    """比较两种方法的结果"""
    print("=" * 70)
    print("橙皮苷逆合成分析结果比较")
    print("=" * 70)

    # 加载结果
    baseline_routes = load_routes(baseline_file)
    flavonoid_routes = load_routes(flavonoid_file)

    print(f"\nBaseline（仅USPTO）结果:")
    print(f"  - 总路线数: {len(baseline_routes)}")

    print(f"\nFlavonoid（USPTO + 黄酮类模板）结果:")
    print(f"  - 总路线数: {len(flavonoid_routes)}")

    # 分析每条路线
    baseline_analysis = [analyze_route(r) for r in baseline_routes]
    flavonoid_analysis = [analyze_route(r) for r in flavonoid_routes]

    # 统计
    baseline_solved = sum(1 for a in baseline_analysis if a['is_solved'])
    flavonoid_solved = sum(1 for a in flavonoid_analysis if a['is_solved'])

    baseline_steps = [a['n_steps'] for a in baseline_analysis if a['is_solved']]
    flavonoid_steps = [a['n_steps'] for a in flavonoid_analysis if a['is_solved']]

    print("\n" + "-" * 70)
    print("解决情况统计:")
    print("-" * 70)
    print(f"Baseline: {baseline_solved}/{len(baseline_routes)} 条路线完全解决")
    print(f"Flavonoid: {flavonoid_solved}/{len(flavonoid_routes)} 条路线完全解决")

    print("\n" + "-" * 70)
    print("解决路线的步骤数统计:")
    print("-" * 70)

    if baseline_steps:
        print(f"Baseline - 最少步骤: {min(baseline_steps)}, 最多步骤: {max(baseline_steps)}, 平均: {sum(baseline_steps)/len(baseline_steps):.1f}")
    else:
        print("Baseline - 无解决路线")

    if flavonoid_steps:
        print(f"Flavonoid - 最少步骤: {min(flavonoid_steps)}, 最多步骤: {max(flavonoid_steps)}, 平均: {sum(flavonoid_steps)/len(flavonoid_steps):.1f}")
    else:
        print("Flavonoid - 无解决路线")

    # 打印详细路线信息
    print("\n" + "=" * 70)
    print("Baseline 路线详情:")
    print("=" * 70)
    for i, (route, analysis) in enumerate(zip(baseline_routes, baseline_analysis)):
        status = "✓ 已解决" if analysis['is_solved'] else "✗ 未解决"
        print(f"路线 {i+1}: {status}, {analysis['n_steps']} 步")

    print("\n" + "=" * 70)
    print("Flavonoid 路线详情:")
    print("=" * 70)
    for i, (route, analysis) in enumerate(zip(flavonoid_routes, flavonoid_analysis)):
        status = "✓ 已解决" if analysis['is_solved'] else "✗ 未解决"
        print(f"路线 {i+1}: {status}, {analysis['n_steps']} 步")

    # 结论
    print("\n" + "=" * 70)
    print("结论:")
    print("=" * 70)

    if flavonoid_solved > baseline_solved:
        improvement = flavonoid_solved - baseline_solved
        print(f"✓ Flavonoid配置比Baseline多找到 {improvement} 条解决路线")
    elif flavonoid_solved == baseline_solved:
        print("两种配置找到的解决路线数相同")
    else:
        print("Baseline配置找到更多解决路线")

    if flavonoid_steps and baseline_steps:
        if min(flavonoid_steps) < min(baseline_steps):
            print(f"✓ Flavonoid配置找到更短的合成路线（{min(flavonoid_steps)} 步 vs {min(baseline_steps)} 步）")
        elif min(flavonoid_steps) == min(baseline_steps):
            print("两种配置的最短路线步骤数相同")
        else:
            print(f"Baseline配置的最短路线更短（{min(baseline_steps)} 步 vs {min(flavonoid_steps)} 步）")

    return {
        'baseline': {
            'total': len(baseline_routes),
            'solved': baseline_solved,
            'min_steps': min(baseline_steps) if baseline_steps else None,
        },
        'flavonoid': {
            'total': len(flavonoid_routes),
            'solved': flavonoid_solved,
            'min_steps': min(flavonoid_steps) if flavonoid_steps else None,
        }
    }


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python compare_results.py <baseline.json> <flavonoid.json>")
        sys.exit(1)

    baseline_file = sys.argv[1]
    flavonoid_file = sys.argv[2]

    compare_results(baseline_file, flavonoid_file)
