#!/usr/bin/env python3
"""
比较baseline、flavonoid和优化配置的逆合成结果
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


def compare_results(baseline_file: str, flavonoid_file: str, optimized_file: str = None):
    """比较三种方法的结果"""
    print("=" * 80)
    print("橙皮苷逆合成分析结果比较")
    print("=" * 80)

    # 加载结果
    baseline_routes = load_routes(baseline_file)
    flavonoid_routes = load_routes(flavonoid_file)
    optimized_routes = load_routes(optimized_file) if optimized_file else None

    print(f"\n配置概览:")
    print(f"  - Baseline（仅USPTO）: C=1.4, 迭代=500, 时间=600秒, 最大变换=10")
    print(f"  - Flavonoid（USPTO + 黄酮类模板）: C=1.4, 迭代=500, 时间=600秒, 最大变换=10")
    if optimized_routes:
        print(f"  - Optimized（优化配置）: C=2.0, 迭代=2000, 时间=1200秒, 最大变换=15")
        print(f"    奖励函数: number of pre-cursors in stock (70%) + fraction in stock (30%)")

    print(f"\n路线总数:")
    print(f"  - Baseline: {len(baseline_routes)}")
    print(f"  - Flavonoid: {len(flavonoid_routes)}")
    if optimized_routes:
        print(f"  - Optimized: {len(optimized_routes)}")

    # 分析每条路线
    baseline_analysis = [analyze_route(r) for r in baseline_routes]
    flavonoid_analysis = [analyze_route(r) for r in flavonoid_routes]
    optimized_analysis = [analyze_route(r) for r in optimized_routes] if optimized_routes else None

    # 统计
    baseline_solved = sum(1 for a in baseline_analysis if a['is_solved'])
    flavonoid_solved = sum(1 for a in flavonoid_analysis if a['is_solved'])
    optimized_solved = sum(1 for a in optimized_analysis if a['is_solved']) if optimized_analysis else None

    baseline_steps = [a['n_steps'] for a in baseline_analysis if a['is_solved']]
    flavonoid_steps = [a['n_steps'] for a in flavonoid_analysis if a['is_solved']]
    optimized_steps = [a['n_steps'] for a in optimized_analysis if a['is_solved']] if optimized_analysis else None

    print("\n" + "-" * 80)
    print("解决情况统计:")
    print("-" * 80)
    print(f"Baseline:   {baseline_solved}/{len(baseline_routes)} 条路线完全解决")
    print(f"Flavonoid:  {flavonoid_solved}/{len(flavonoid_routes)} 条路线完全解决")
    if optimized_routes:
        print(f"Optimized:  {optimized_solved}/{len(optimized_routes)} 条路线完全解决")

    print("\n" + "-" * 80)
    print("解决路线的步骤数统计:")
    print("-" * 80)

    if baseline_steps:
        print(f"Baseline  - 最少步骤: {min(baseline_steps)}, 最多步骤: {max(baseline_steps)}, 平均: {sum(baseline_steps)/len(baseline_steps):.1f}")
    else:
        print("Baseline  - 无解决路线")

    if flavonoid_steps:
        print(f"Flavonoid - 最少步骤: {min(flavonoid_steps)}, 最多步骤: {max(flavonoid_steps)}, 平均: {sum(flavonoid_steps)/len(flavonoid_steps):.1f}")
    else:
        print("Flavonoid - 无解决路线")

    if optimized_routes:
        if optimized_steps:
            print(f"Optimized - 最少步骤: {min(optimized_steps)}, 最多步骤: {max(optimized_steps)}, 平均: {sum(optimized_steps)/len(optimized_steps):.1f}")
        else:
            print("Optimized - 无解决路线")

    # 打印详细路线信息
    print("\n" + "=" * 80)
    print("Baseline 路线详情:")
    print("=" * 80)
    for i, (route, analysis) in enumerate(zip(baseline_routes, baseline_analysis)):
        status = "✓ 已解决" if analysis['is_solved'] else "✗ 未解决"
        print(f"路线 {i+1}: {status}, {analysis['n_steps']} 步")

    print("\n" + "=" * 80)
    print("Flavonoid 路线详情:")
    print("=" * 80)
    for i, (route, analysis) in enumerate(zip(flavonoid_routes, flavonoid_analysis)):
        status = "✓ 已解决" if analysis['is_solved'] else "✗ 未解决"
        print(f"路线 {i+1}: {status}, {analysis['n_steps']} 步")

    if optimized_routes:
        print("\n" + "=" * 80)
        print("Optimized 路线详情:")
        print("=" * 80)
        for i, (route, analysis) in enumerate(zip(optimized_routes, optimized_analysis)):
            status = "✓ 已解决" if analysis['is_solved'] else "✗ 未解决"
            print(f"路线 {i+1}: {status}, {analysis['n_steps']} 步")

    # 结论
    print("\n" + "=" * 80)
    print("结论:")
    print("=" * 80)

    results = {
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

    if optimized_routes:
        results['optimized'] = {
            'total': len(optimized_routes),
            'solved': optimized_solved,
            'min_steps': min(optimized_steps) if optimized_steps else None,
        }

    # 比较Flavonoid和Baseline
    if flavonoid_solved > baseline_solved:
        improvement = flavonoid_solved - baseline_solved
        print(f"✓ Flavonoid配置比Baseline多找到 {improvement} 条解决路线")
    elif flavonoid_solved == baseline_solved:
        print("Baseline和Flavonoid配置找到的解决路线数相同")
    else:
        print("Baseline配置比Flavonoid找到更多解决路线")

    # 比较Optimized和其他配置
    if optimized_routes:
        if optimized_solved > max(baseline_solved, flavonoid_solved):
            print(f"✓ 优化配置表现最佳，找到 {optimized_solved} 条解决路线")
            if optimized_solved > flavonoid_solved:
                print(f"  比Flavonoid多 {optimized_solved - flavonoid_solved} 条解决路线")
            if optimized_solved > baseline_solved:
                print(f"  比Baseline多 {optimized_solved - baseline_solved} 条解决路线")
        elif optimized_solved == max(baseline_solved, flavonoid_solved):
            print("优化配置与其他最佳配置找到相同数量的解决路线")
        else:
            print(f"优化配置找到 {optimized_solved} 条解决路线，少于其他配置")

    # 比较最短路线
    all_configs = [
        ('Baseline', baseline_steps),
        ('Flavonoid', flavonoid_steps),
    ]
    if optimized_routes:
        all_configs.append(('Optimized', optimized_steps))

    configs_with_solutions = [(name, steps) for name, steps in all_configs if steps]

    if configs_with_solutions:
        print("\n最短路线比较:")
        for name, steps in configs_with_solutions:
            print(f"  {name}: {min(steps)} 步")

        best_config = min(configs_with_solutions, key=lambda x: min(x[1]))
        if len(configs_with_solutions) > 1:
            print(f"✓ {best_config[0]} 配置找到最短的合成路线")

    return results


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("用法: python compare_all_results.py <baseline.json> <flavonoid.json> [optimized.json]")
        sys.exit(1)

    baseline_file = sys.argv[1]
    flavonoid_file = sys.argv[2]
    optimized_file = sys.argv[3] if len(sys.argv) == 4 else None

    compare_results(baseline_file, flavonoid_file, optimized_file)
