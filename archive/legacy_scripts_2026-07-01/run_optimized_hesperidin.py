#!/usr/bin/env python3
"""
使用优化配置对橙皮苷进行逆合成分析
优化内容:
1. 使用多重奖励函数优先考虑stock可用性
2. 增加迭代次数和时间限制
3. 提高探索参数C以发现更长的解决路线
4. 增加最大变换步数以允许更长的合成路径
"""
import json
import sys
import time
from pathlib import Path
from datetime import datetime

from aizynthfinder.aizynthfinder import AiZynthFinder

# 橙皮苷 SMILES
HESPERIDIN_SMILES = "COc1ccc(C2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)[C@H](O)[C@@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O"

CONFIG_FILE = "/home/ljx/retro_synthesis/config/hesperidin_optimized.yml"
OUTPUT_DIR = Path("/home/ljx/retro_synthesis/outputs")


def run_retrosynthesis():
    """运行逆合成分析"""
    print("=" * 60)
    print("橙皮苷逆合成分析 - 优化配置")
    print("=" * 60)
    print(f"目标分子 SMILES: {HESPERIDIN_SMILES}")
    print(f"配置文件: {CONFIG_FILE}")
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    print("优化参数:")
    print("  - C=2.0 (提高探索)")
    print("  - 迭代次数=2000 (原500)")
    print("  - 时间限制=1200秒 (原600秒)")
    print("  - 最大变换=15 (原10)")
    print("  - 奖励函数: number of pre-cursors in stock (70%) + fraction in stock (30%)")
    print()

    # 初始化AiZynthFinder
    print("正在初始化AiZynthFinder...")
    finder = AiZynthFinder(configfile=CONFIG_FILE)

    # 选择expansion policy
    print(f"可用的expansion policies: {finder.expansion_policy.items}")
    finder.expansion_policy.select(["uspto", "ringbreaker", "flavonoid", "flavonoid_structural", "flavonoid_biosynthesis"])
    print(f"已选择的expansion policies: {finder.expansion_policy.selection}")

    # 选择库存数据库
    print(f"可用的库存: {finder.stock.items}")
    finder.stock.select(["zinc", "flavonoid", "flavonoid_structural", "flavonoid_biosynthesis"])
    print(f"已选择的库存: {finder.stock.selection}")

    # 设置目标分子
    print("正在设置目标分子...")
    finder.target_smiles = HESPERIDIN_SMILES

    # 准备搜索树
    print("正在准备搜索树...")
    start_time = time.time()
    finder.prepare_tree()

    # 执行搜索
    print("正在执行MCTS搜索...")
    finder.tree_search()

    # 构建路线
    print("正在构建合成路线...")
    finder.build_routes()

    elapsed = time.time() - start_time

    # 提取统计信息
    stats = finder.extract_statistics()

    print(f"\n搜索完成! 耗时: {elapsed:.2f}秒")
    print(f"找到路线数: {len(finder.routes)}")
    print(f"统计信息: {json.dumps(stats, indent=2)}")

    # 保存结果
    timestamp = datetime.now().strftime("%y-%m-%d")
    output_file = OUTPUT_DIR / f"{timestamp}-hesperidin_optimized.json"

    routes_data = []
    for i, route in enumerate(finder.routes):
        reaction_tree = route['reaction_tree']
        route_dict = reaction_tree.to_dict()
        route_dict['score'] = route.get('score', {})
        route_dict['is_solved'] = route.get('route_metadata', {}).get('is_solved', False)
        routes_data.append(route_dict)

    with open(output_file, 'w') as f:
        json.dump(routes_data, f, indent=2)

    print(f"\n结果已保存到: {output_file}")
    print(f"共保存 {len(routes_data)} 条路线")

    # 打印路线摘要
    print("\n" + "=" * 60)
    print("路线摘要:")
    print("=" * 60)
    for i, route in enumerate(finder.routes):
        is_solved = route.get('route_metadata', {}).get('is_solved', False)
        reaction_tree = route['reaction_tree']
        n_steps = len(list(reaction_tree.reactions()))
        print(f"路线 {i+1}: 已解决={is_solved}, 步骤数={n_steps}")

    return output_file, stats


if __name__ == "__main__":
    output_file, stats = run_retrosynthesis()
    print(f"\n完成! 输出文件: {output_file}")
