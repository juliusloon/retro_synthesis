#!/usr/bin/env python3
"""
将AiZynthFinder的JSON输出转换为可读的Markdown摘要报告

用法:
    python json2md.py <input.json> [output.md]

如果不指定输出文件，将自动生成同名的.md文件
"""
import json
import sys
from pathlib import Path
from typing import List, Dict, Any, Optional


def count_reactions(node: Dict[str, Any]) -> int:
    """计算路线中的反应步骤数"""
    count = 0
    if isinstance(node, dict):
        if node.get('type') == 'reaction':
            count += 1
        for child in node.get('children', []):
            count += count_reactions(child)
    return count


def get_leaf_molecules(node: Dict[str, Any]) -> List[Dict[str, str]]:
    """获取路线中的叶子节点（起始原料）"""
    leaves = []
    if isinstance(node, dict):
        children = node.get('children', [])
        if node.get('type') == 'mol' and (not children or len(children) == 0):
            leaves.append({
                'smiles': node.get('smiles', ''),
                'in_stock': node.get('in_stock', False)
            })
        else:
            for child in children:
                leaves.extend(get_leaf_molecules(child))
    return leaves


def get_reaction_types(node: Dict[str, Any]) -> List[str]:
    """获取路线中使用的反应类型（分类）"""
    types = []
    if isinstance(node, dict):
        if node.get('type') == 'reaction':
            classification = node.get('metadata', {}).get('classification', 'Unknown')
            policy = node.get('metadata', {}).get('policy_name', 'Unknown')
            types.append(f"{classification} ({policy})")
        for child in node.get('children', []):
            types.extend(get_reaction_types(child))
    return types


def get_policies_used(node: Dict[str, Any]) -> Dict[str, int]:
    """统计使用的策略及次数"""
    policies = {}
    if isinstance(node, dict):
        if node.get('type') == 'reaction':
            policy = node.get('metadata', {}).get('policy_name', 'Unknown')
            policies[policy] = policies.get(policy, 0) + 1
        for child in node.get('children', []):
            child_policies = get_policies_used(child)
            for k, v in child_policies.items():
                policies[k] = policies.get(k, 0) + v
    return policies


def is_route_solved(node: Dict[str, Any]) -> bool:
    """检查路线是否完全解决（所有起始原料都在库存中）"""
    leaves = get_leaf_molecules(node)
    if not leaves:
        return False
    return all(leaf['in_stock'] for leaf in leaves)


def analyze_route(route: Dict[str, Any], index: int) -> Dict[str, Any]:
    """分析单条路线，返回摘要信息"""
    n_steps = count_reactions(route)
    leaves = get_leaf_molecules(route)
    reaction_types = get_reaction_types(route)
    policies = get_policies_used(route)

    solved = route.get('is_solved', False)
    if not solved:
        solved = is_route_solved(route)

    score = route.get('score', {})

    in_stock_count = sum(1 for leaf in leaves if leaf['in_stock'])
    total_leaves = len(leaves)

    return {
        'index': index,
        'n_steps': n_steps,
        'is_solved': solved,
        'score': score,
        'leaves': leaves,
        'in_stock_count': in_stock_count,
        'total_leaves': total_leaves,
        'stock_fraction': in_stock_count / total_leaves if total_leaves > 0 else 0,
        'reaction_types': reaction_types,
        'policies': policies,
        'target_smiles': route.get('smiles', '')
    }


def generate_markdown(routes: List[Dict[str, Any]], input_file: str) -> str:
    """生成Markdown报告"""
    lines = []

    # 标题
    lines.append("# 逆合成路线分析报告")
    lines.append("")
    lines.append(f"**输入文件**: `{input_file}`")
    lines.append(f"**路线总数**: {len(routes)}")
    lines.append("")

    # 分析所有路线
    analyses = [analyze_route(route, i+1) for i, route in enumerate(routes)]

    # 总体统计
    solved_count = sum(1 for a in analyses if a['is_solved'])
    step_counts = [a['n_steps'] for a in analyses]

    lines.append("## 总体统计")
    lines.append("")
    lines.append(f"| 指标 | 值 |")
    lines.append(f"|------|-----|")
    lines.append(f"| 完全解决路线数 | {solved_count}/{len(routes)} |")
    lines.append(f"| 最短路线 | {min(step_counts)} 步 |")
    lines.append(f"| 最长路线 | {max(step_counts)} 步 |")
    lines.append(f"| 平均路线长度 | {sum(step_counts)/len(step_counts):.1f} 步 |")
    lines.append("")

    # 分组统计
    short_routes = [a for a in analyses if a['n_steps'] <= 5]
    medium_routes = [a for a in analyses if 5 < a['n_steps'] <= 15]
    long_routes = [a for a in analyses if a['n_steps'] > 15]

    lines.append("### 路线长度分布")
    lines.append("")
    lines.append(f"| 类别 | 数量 | 占比 |")
    lines.append(f"|------|------|------|")
    lines.append(f"| 短路线 (≤5步) | {len(short_routes)} | {len(short_routes)/len(analyses)*100:.1f}% |")
    lines.append(f"| 中等路线 (6-15步) | {len(medium_routes)} | {len(medium_routes)/len(analyses)*100:.1f}% |")
    lines.append(f"| 长路线 (>15步) | {len(long_routes)} | {len(long_routes)/len(analyses)*100:.1f}% |")
    lines.append("")

    # 路线详情表格
    lines.append("## 路线详情")
    lines.append("")
    lines.append("| # | 步骤数 | 已解决 | 库存分数 | 评分 | 起始原料数量 | 在库数量 |")
    lines.append("|---|--------|--------|----------|------|--------------|----------|")

    for a in analyses:
        solved_icon = "✅" if a['is_solved'] else "❌"
        score_str = ", ".join([f"{k}: {v}" for k, v in a['score'].items()]) if a['score'] else "N/A"
        lines.append(f"| {a['index']} | {a['n_steps']} | {solved_icon} | {a['stock_fraction']:.2f} | {score_str} | {a['total_leaves']} | {a['in_stock_count']} |")

    lines.append("")

    # 详细路线信息
    lines.append("## 路线详细信息")
    lines.append("")

    for a in analyses:
        lines.append(f"### 路线 {a['index']}")
        lines.append("")
        lines.append(f"- **步骤数**: {a['n_steps']}")
        lines.append(f"- **是否解决**: {'✅ 是' if a['is_solved'] else '❌ 否'}")
        lines.append(f"- **起始原料**: {a['total_leaves']} 个 (在库: {a['in_stock_count']})")
        lines.append(f"- **评分**: {a['score']}")
        lines.append("")

        # 使用的策略
        if a['policies']:
            lines.append("**使用的策略**:")
            for policy, count in sorted(a['policies'].items(), key=lambda x: -x[1]):
                lines.append(f"- {policy}: {count} 次")
            lines.append("")

        # 起始原料列表
        if a['leaves']:
            lines.append("**起始原料**:")
            for leaf in a['leaves']:
                stock_icon = "✅" if leaf['in_stock'] else "❌"
                lines.append(f"- {stock_icon} `{leaf['smiles']}`")
            lines.append("")

        lines.append("---")
        lines.append("")

    return "\n".join(lines)


def main():
    if len(sys.argv) < 2:
        print("用法: python json2md.py <input.json> [output.md]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.replace('.json', '.md')

    # 加载JSON
    with open(input_file, 'r') as f:
        routes = json.load(f)

    print(f"正在分析 {input_file}...")
    print(f"找到 {len(routes)} 条路线")

    # 生成Markdown
    md_content = generate_markdown(routes, input_file)

    # 保存文件
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(md_content)

    print(f"报告已保存到: {output_file}")


if __name__ == "__main__":
    main()
