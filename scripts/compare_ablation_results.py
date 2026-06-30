#!/usr/bin/env python3
"""
比较消融实验结果，输出增强评价指标。

增强指标:
  - AiZynth solved（AiZynthFinder 原生 is_solved）
  - map-valid solved（所有反应步骤通过 atom-map 守恒审计）
  - strict-stock + map-valid solved（真实库存 + map 守恒）
  - 每条路线的 stock 来源分布
  - 自定义模板贡献（custom step 数）
  - 路线化学合理性（是否先断糖苷键）
  - stock 层级分析（strict_buyable, trusted_intermediate, virtual_bridge）
  - 自定义模板家族统计
  - 关键中间体到达分析（aglycone, chalcone, sugar donor）

用法:
    python scripts/compare_ablation_results.py
"""
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import inchi

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "ablation"
STOCK_LAYERS_DIR = PROJECT_ROOT / "templates" / "stock_layers"

# Map retention 阈值：目标侧 map numbers 在前体侧保留的比例
MAP_RETENTION_THRESHOLD = 0.8


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


def load_reaction_families() -> dict:
    """加载反应家族配置，返回分类到家族名称的映射。"""
    families_file = PROJECT_ROOT / "config" / "reaction_families.json"
    if not families_file.exists():
        return {}
    with open(families_file) as f:
        families = json.load(f)
    return {family["classification"]: family["name"] for family in families}


def load_stock_metadata() -> dict:
    """加载 stock 元数据，返回 InChIKey 到 stock 层级的映射。"""
    metadata_file = STOCK_LAYERS_DIR / "stock_layers_metadata.csv"
    if not metadata_file.exists():
        return {}
    inchikey_to_layer = {}
    with open(metadata_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            inchikey_to_layer[row["inchikey"]] = row["stock_layer"]
    return inchikey_to_layer


def load_routes(filepath: str) -> list:
    with open(filepath) as f:
        return json.load(f)


def smiles_to_inchikey(smiles: str) -> str:
    """从 SMILES 计算 InChIKey。"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        return inchi.MolToInchiKey(mol)
    except Exception:
        return ""


def compute_route_map_retention(route_dict: dict) -> float:
    """计算整条路线的最小 atom-map retention ratio。

    AiZynthFinder stores retro reactions as:
      current target molecule >> retrosynthetic precursors

    For each reaction node, compute:
      target_map_numbers intersect precursor_map_numbers / target_map_numbers
    返回所有 reaction 节点中的最小值。
    若无 reaction 节点，返回 1.0。
    """
    min_ratio = 1.0

    def walk(node):
        nonlocal min_ratio
        if not isinstance(node, dict):
            return
        if node.get('type') == 'reaction':
            meta = node.get('metadata', {})
            mapped = meta.get('mapped_reaction_smiles', '')
            if '>>' in mapped:
                tgt_side, prec_side = mapped.split('>>', 1)
                tgt_maps = set(re.findall(r':(\d+)\]', tgt_side))
                prec_maps = set(re.findall(r':(\d+)\]', prec_side))
                if tgt_maps:
                    ratio = len(tgt_maps & prec_maps) / len(tgt_maps)
                    min_ratio = min(min_ratio, ratio)
        for child in node.get('children', []):
            walk(child)

    walk(route_dict)
    return min_ratio


def analyze_route(route_dict: dict, reaction_families: dict = None) -> dict:
    """分析单条路线的详细信息。"""
    if reaction_families is None:
        reaction_families = {}

    def walk_tree(node, depth=0):
        """递归遍历路线树，提取步骤信息。"""
        steps = []
        if isinstance(node, dict):
            if node.get('type') == 'reaction':
                metadata = node.get('metadata', {})
                step_info = {
                    'policy_name': metadata.get('policy_name', ''),
                    'template_code': metadata.get('template_code', ''),
                    'classification': metadata.get('classification', ''),
                    'template': metadata.get('template', ''),
                    'depth': depth,
                }
                steps.append(step_info)
            for child in node.get('children', []):
                steps.extend(walk_tree(child, depth + 1))
        return steps

    def count_reactions(node):
        """只统计 type=='reaction' 的节点，不把根分子算作一步。"""
        count = 0
        if isinstance(node, dict):
            if node.get('type') == 'reaction':
                count += 1
            for child in node.get('children', []):
                count += count_reactions(child)
        return count

    def is_solved_recursive(node):
        """递归检查：所有叶子 mol 节点必须满足 in_stock == True。"""
        if isinstance(node, dict):
            children = node.get('children', [])
            if not children:
                if node.get('type') == 'mol':
                    return node.get('in_stock', False) is True
                return False
            return all(is_solved_recursive(child) for child in children)
        return False

    def collect_leaf_smiles(node):
        """收集叶子节点的 SMILES（起始原料）。"""
        leaves = []
        if isinstance(node, dict):
            if node.get('type') == 'mol':
                children = node.get('children', [])
                if not children or len(children) == 0:
                    leaves.append(node.get('smiles', ''))
                else:
                    for child in children:
                        leaves.extend(collect_leaf_smiles(child))
            else:
                for child in node.get('children', []):
                    leaves.extend(collect_leaf_smiles(child))
        return leaves

    def collect_leaf_inchikeys(node):
        """收集叶子节点的 InChIKey。"""
        leaves = []
        if isinstance(node, dict):
            if node.get('type') == 'mol':
                children = node.get('children', [])
                if not children or len(children) == 0:
                    inchikey = node.get('inchikey', '')
                    if not inchikey:
                        inchikey = smiles_to_inchikey(node.get('smiles', ''))
                    if inchikey:
                        leaves.append(inchikey)
                else:
                    for child in children:
                        leaves.extend(collect_leaf_inchikeys(child))
            else:
                for child in node.get('children', []):
                    leaves.extend(collect_leaf_inchikeys(child))
        return leaves

    def has_glycoside_bond_cleavage(steps):
        """检查路线中是否有糖苷键断裂步骤。"""
        glyco_keywords = ['glycosylation', 'o_glycosylation', 'c_glycosylation',
                          'glycoside', 'sugar', 'glucoside']
        for step in steps:
            cls = step.get('classification', '').lower()
            template = step.get('template', '').lower()
            if any(kw in cls or kw in template for kw in glyco_keywords):
                return True
        return False

    def has_chalcone_retro_aldol(steps):
        """检查路线中是否有 chalcone retro aldol 步骤。"""
        chalcone_keywords = ['chalcone_retro_aldol', 'chalcone_aldol', 'retro_aldol']
        for step in steps:
            cls = step.get('classification', '').lower()
            template = step.get('template', '').lower()
            if any(kw in cls or kw in template for kw in chalcone_keywords):
                return True
        return False

    def analyze_custom_family_counts(steps):
        """统计自定义模板家族的使用次数。"""
        family_counts = defaultdict(int)
        for step in steps:
            pname = step.get('policy_name', 'unknown')
            if pname not in ('uspto', 'ringbreaker'):
                cls = step.get('classification', '')
                family_name = reaction_families.get(cls, cls)
                family_counts[family_name] += 1
        return dict(family_counts)

    def get_first_custom_step_family(steps):
        """获取第一个自定义模板步骤的家族。"""
        for step in steps:
            pname = step.get('policy_name', 'unknown')
            if pname not in ('uspto', 'ringbreaker'):
                cls = step.get('classification', '')
                return reaction_families.get(cls, cls)
        return None

    def check_route_reaches_intermediates(steps):
        """检查路线是否到达了关键中间体。"""
        # 检查是否到达 aglycone（通过 glycoside cleavage 或相关分类）
        glyco_keywords = ['glycosylation', 'o_glycosylation', 'c_glycosylation',
                          'glycoside', 'sugar', 'glucoside']
        chalcone_keywords = ['chalcone_retro_aldol', 'chalcone_aldol', 'retro_aldol',
                            'flavanone_chalcone_opening']
        sugar_keywords = ['glycosylation', 'o_glycosylation', 'sugar']

        reaches_aglycone = False
        reaches_chalcone = False
        reaches_sugar_donor = False

        for step in steps:
            cls = step.get('classification', '').lower()
            template = step.get('template', '').lower()

            if any(kw in cls or kw in template for kw in glyco_keywords):
                reaches_aglycone = True
            if any(kw in cls or kw in template for kw in chalcone_keywords):
                reaches_chalcone = True
            if any(kw in cls or kw in template for kw in sugar_keywords):
                reaches_sugar_donor = True

        return reaches_aglycone, reaches_chalcone, reaches_sugar_donor

    n_steps = count_reactions(route_dict)
    # 优先使用 AiZynthFinder 写入的顶层 is_solved 字段
    aizynth_solved = route_dict.get('is_solved')
    if aizynth_solved is None:
        aizynth_solved = is_solved_recursive(route_dict)
    # 递归检查（用于交叉验证）
    recursive_solved = is_solved_recursive(route_dict)
    steps = walk_tree(route_dict)
    leaf_smiles = collect_leaf_smiles(route_dict)
    leaf_inchikeys = collect_leaf_inchikeys(route_dict)

    # Map retention 审计
    map_retention = compute_route_map_retention(route_dict)
    map_valid = map_retention >= MAP_RETENTION_THRESHOLD

    # Policy 来源统计
    policy_counts = defaultdict(int)
    custom_steps = 0
    for step in steps:
        pname = step.get('policy_name', 'unknown')
        policy_counts[pname] += 1
        if pname not in ('uspto', 'ringbreaker'):
            custom_steps += 1

    # 自定义模板家族统计
    custom_family_counts = analyze_custom_family_counts(steps)
    first_custom_step_family = get_first_custom_step_family(steps)

    # 关键中间体到达分析
    reaches_aglycone, reaches_chalcone, reaches_sugar_donor = check_route_reaches_intermediates(steps)

    return {
        'n_steps': n_steps,
        'aizynth_solved': aizynth_solved,
        'recursive_solved': recursive_solved,
        'map_retention': map_retention,
        'map_valid': map_valid,
        'is_solved': aizynth_solved and map_valid,  # 真正有效：AiZynth solved + map 守恒
        'steps': steps,
        'leaf_smiles': leaf_smiles,
        'leaf_inchikeys': leaf_inchikeys,
        'policy_counts': dict(policy_counts),
        'custom_steps': custom_steps,
        'has_glycoside_cleavage': has_glycoside_bond_cleavage(steps),
        'has_chalcone_retro_aldol': has_chalcone_retro_aldol(steps),
        'custom_family_counts': custom_family_counts,
        'first_custom_step_family': first_custom_step_family,
        'reaches_aglycone': reaches_aglycone,
        'reaches_chalcone': reaches_chalcone,
        'reaches_sugar_donor': reaches_sugar_donor,
    }


def analyze_stock_layers(leaf_inchikeys: list, stock_layers: dict) -> dict:
    """分析叶子节点的 stock 层级分布。"""
    layer_counts = {layer: 0 for layer in stock_layers.keys()}
    layer_counts["unknown"] = 0

    for inchikey in leaf_inchikeys:
        found = False
        for layer_name, inchikeys in stock_layers.items():
            if inchikey in inchikeys:
                layer_counts[layer_name] += 1
                found = True
                break
        if not found:
            layer_counts["unknown"] += 1

    return layer_counts


def route_uses_stock_layer(leaf_inchikeys: list, stock_layers: dict, layer_name: str) -> bool:
    """检查路线叶子是否使用某个 stock 层级。"""
    layer_keys = stock_layers.get(layer_name, set())
    return any(inchikey in layer_keys for inchikey in leaf_inchikeys)


def compare_ablation():
    """比较所有消融实验结果。"""
    # 加载 stock 层级信息和反应家族
    stock_layers = load_stock_layer_inchikeys()
    reaction_families = load_reaction_families()
    stock_metadata = load_stock_metadata()

    result_files = sorted(OUTPUT_DIR.glob("*.json"))
    result_files = [f for f in result_files if f.name not in (
        "ablation_summary.json", "ablation_report.json", "route_gap_worklist.json")]

    if not result_files:
        print("错误: 未找到消融实验结果文件")
        print(f"请先运行: python scripts/run_ablation_experiments.py")
        sys.exit(1)

    print("=" * 90)
    print("橙皮苷逆合成消融实验结果比较")
    print("=" * 90)

    all_results = {}

    for filepath in result_files:
        name = filepath.stem
        routes = load_routes(str(filepath))
        analyses = [analyze_route(r, reaction_families) for r in routes]

        n_total = len(routes)
        n_aizynth_solved = sum(1 for a in analyses if a['aizynth_solved'])
        n_map_valid = sum(1 for a in analyses if a['aizynth_solved'] and a['map_valid'])
        n_effective = sum(1 for a in analyses if a['is_solved'])
        solved_analyses = [a for a in analyses if a['is_solved']]
        bridge_closed_analyses = [
            a for a in solved_analyses
            if route_uses_stock_layer(a['leaf_inchikeys'], stock_layers, "virtual_bridge")
        ]
        non_virtual_solved_analyses = [
            a for a in solved_analyses
            if not route_uses_stock_layer(a['leaf_inchikeys'], stock_layers, "virtual_bridge")
        ]
        n_bridge_closed = len(bridge_closed_analyses)
        n_non_virtual_effective = len(non_virtual_solved_analyses)

        # 步骤数统计（仅对有效 solved 路线）
        step_counts = [a['n_steps'] for a in solved_analyses]

        # 自定义模板贡献（solved + all routes）
        custom_step_total_solved = sum(a['custom_steps'] for a in solved_analyses)
        custom_step_total_all = sum(a['custom_steps'] for a in analyses)

        # 糖苷键断裂（solved + all routes）
        glyco_count_solved = sum(1 for a in solved_analyses if a['has_glycoside_cleavage'])
        glyco_count_all = sum(1 for a in analyses if a['has_glycoside_cleavage'])

        # Chalcone retro aldol（solved + all routes）
        chalcone_retro_aldol_count_solved = sum(1 for a in solved_analyses if a['has_chalcone_retro_aldol'])
        chalcone_retro_aldol_count_all = sum(1 for a in analyses if a['has_chalcone_retro_aldol'])

        # Policy 来源（all routes）
        all_policies = defaultdict(int)
        for a in analyses:
            for pname, cnt in a['policy_counts'].items():
                all_policies[pname] += cnt

        # Map retention 分布
        map_retentions = [a['map_retention'] for a in analyses]

        # Stock 层级分析（solved routes）
        stock_layer_solved = {layer: 0 for layer in stock_layers.keys()}
        stock_layer_solved["unknown"] = 0
        for a in solved_analyses:
            layer_counts = analyze_stock_layers(a['leaf_inchikeys'], stock_layers)
            for layer, count in layer_counts.items():
                if count > 0:
                    stock_layer_solved[layer] += 1

        # 自定义模板家族统计（solved + all routes）
        solved_custom_family_counts = defaultdict(int)
        for a in solved_analyses:
            for family, count in a['custom_family_counts'].items():
                solved_custom_family_counts[family] += count
        all_custom_family_counts = defaultdict(int)
        for a in analyses:
            for family, count in a['custom_family_counts'].items():
                all_custom_family_counts[family] += count

        # 关键中间体到达统计（solved + all routes）
        reaches_aglycone_count_solved = sum(1 for a in solved_analyses if a['reaches_aglycone'])
        reaches_chalcone_count_solved = sum(1 for a in solved_analyses if a['reaches_chalcone'])
        reaches_sugar_donor_count_solved = sum(1 for a in solved_analyses if a['reaches_sugar_donor'])
        reaches_aglycone_count_all = sum(1 for a in analyses if a['reaches_aglycone'])
        reaches_chalcone_count_all = sum(1 for a in analyses if a['reaches_chalcone'])
        reaches_sugar_donor_count_all = sum(1 for a in analyses if a['reaches_sugar_donor'])

        all_results[name] = {
            'n_total': n_total,
            'n_aizynth_solved': n_aizynth_solved,
            'n_map_valid': n_map_valid,
            'n_effective': n_effective,
            'n_bridge_closed': n_bridge_closed,
            'n_non_virtual_effective': n_non_virtual_effective,
            'step_counts': step_counts,
            'custom_step_total': custom_step_total_solved,
            'custom_step_total_all': custom_step_total_all,
            'glyco_count': glyco_count_solved,
            'glyco_count_all': glyco_count_all,
            'chalcone_retro_aldol_count': chalcone_retro_aldol_count_solved,
            'chalcone_retro_aldol_count_all': chalcone_retro_aldol_count_all,
            'policy_counts': dict(all_policies),
            'map_retentions': map_retentions,
            'analyses': analyses,
            'stock_layer_solved': stock_layer_solved,
            'custom_family_counts': dict(solved_custom_family_counts),
            'all_custom_family_counts': dict(all_custom_family_counts),
            'reaches_aglycone_count': reaches_aglycone_count_solved,
            'reaches_chalcone_count': reaches_chalcone_count_solved,
            'reaches_sugar_donor_count': reaches_sugar_donor_count_solved,
            'reaches_aglycone_count_all': reaches_aglycone_count_all,
            'reaches_chalcone_count_all': reaches_chalcone_count_all,
            'reaches_sugar_donor_count_all': reaches_sugar_donor_count_all,
        }

    # 三列汇总表
    print(f"\n{'实验':<30} {'AiZynth':>8} {'map-valid':>10} {'有效solved':>10} {'bridge':>8} {'non-virt':>8} {'最短步':>6}")
    print("-" * 88)
    for name, r in all_results.items():
        min_step = min(r['step_counts']) if r['step_counts'] else '-'
        print(f"{name:<30} {r['n_aizynth_solved']:>8} {r['n_map_valid']:>10} "
              f"{r['n_effective']:>10} {r['n_bridge_closed']:>8} "
              f"{r['n_non_virtual_effective']:>8} {str(min_step):>6}")

    # 详细分析
    print(f"\n{'='*90}")
    print("详细分析")
    print(f"{'='*90}")

    for name, r in all_results.items():
        print(f"\n--- {name} ---")
        print(f"  AiZynth solved: {r['n_aizynth_solved']}/{r['n_total']}")
        print(f"  map-valid solved: {r['n_map_valid']}/{r['n_total']}")
        print(f"  有效 solved (AiZynth + map): {r['n_effective']}/{r['n_total']}")
        print(f"  bridge-closed solved (uses virtual_bridge): {r['n_bridge_closed']}/{r['n_total']}")
        print(f"  non-virtual effective solved: {r['n_non_virtual_effective']}/{r['n_total']}")

        # Map retention 分布
        rets = r['map_retentions']
        if rets:
            print(f"  Map retention: min={min(rets):.2f}, max={max(rets):.2f}, "
                  f"avg={sum(rets)/len(rets):.2f}")
            zero_ret = sum(1 for x in rets if x == 0.0)
            if zero_ret > 0:
                print(f"  ⚠ {zero_ret}/{len(rets)} 条路线 map retention = 0（scaffold collapse）")

        if r['step_counts']:
            print(f"  有效路线步骤数: min={min(r['step_counts'])}, "
                  f"max={max(r['step_counts'])}, "
                  f"avg={sum(r['step_counts'])/len(r['step_counts']):.1f}")

        if r['policy_counts']:
            print(f"  Policy 来源: {r['policy_counts']}")

        print(f"  自定义模板步骤总数: solved={r['custom_step_total']} / all={r['custom_step_total_all']}")
        print(f"  含糖苷键断裂的路线: solved={r['glyco_count']} / all={r['glyco_count_all']}")
        print(f"  含 chalcone retro aldol: solved={r['chalcone_retro_aldol_count']} / all={r['chalcone_retro_aldol_count_all']}")

        # Stock 层级分析
        print(f"  Stock 层级 solved 分布:")
        for layer, count in r['stock_layer_solved'].items():
            if count > 0:
                print(f"    {layer}: {count} 条路线")

        # 自定义模板家族统计（all routes）
        if r['all_custom_family_counts']:
            print(f"  自定义模板家族统计 (all routes):")
            for family, count in sorted(r['all_custom_family_counts'].items(), key=lambda x: -x[1]):
                print(f"    {family}: {count} 次")
        if r['custom_family_counts']:
            print(f"  自定义模板家族统计 (solved only):")
            for family, count in sorted(r['custom_family_counts'].items(), key=lambda x: -x[1]):
                print(f"    {family}: {count} 次")

        # 关键中间体到达统计
        print(f"  关键中间体到达 (solved / all):")
        print(f"    aglycone: {r['reaches_aglycone_count']} / {r['reaches_aglycone_count_all']}")
        print(f"    chalcone: {r['reaches_chalcone_count']} / {r['reaches_chalcone_count_all']}")
        print(f"    sugar donor: {r['reaches_sugar_donor_count']} / {r['reaches_sugar_donor_count_all']}")

        # 列出每条路线的详情
        for i, a in enumerate(r['analyses']):
            status = "✓" if a['is_solved'] else "✗"
            aiz = "✓" if a['aizynth_solved'] else "✗"
            mval = "✓" if a['map_valid'] else "✗"
            ret = f"{a['map_retention']:.2f}"
            policies = ', '.join(f"{k}:{v}" for k, v in a['policy_counts'].items())
            custom_families = ', '.join(f"{k}:{v}" for k, v in a['custom_family_counts'].items()) if a['custom_family_counts'] else "none"
            print(f"  路线 {i+1}: [{status}] aizynth={aiz} map={mval}({ret}) "
                  f"{a['n_steps']}步 custom={a['custom_steps']} policies=[{policies}]")
            if a['custom_family_counts']:
                print(f"    自定义家族: [{custom_families}]")
            print(f"    到达: aglycone={a['reaches_aglycone']} chalcone={a['reaches_chalcone']} sugar={a['reaches_sugar_donor']}")

    # 消融对比分析
    print(f"\n{'='*90}")
    print("消融对比分析")
    print(f"{'='*90}")

    def get(name):
        return all_results.get(name, {})

    baseline_zinc = get("baseline_zinc")
    baseline_strict = get("baseline_strict")
    flavonoid_zinc = get("flavonoid_zinc")
    flavonoid_strict = get("flavonoid_strict")
    custom_only = get("custom_only_strict")
    custom_only_virtual = get("custom_only_virtual_bridge")
    flavonoid_virtual = get("flavonoid_virtual_bridge")

    # 1. ZINC vs Strict stock 的影响 (baseline)
    if baseline_zinc and baseline_strict:
        print(f"\n[1] Stock 影响 (Baseline):")
        print(f"  ZINC:      {baseline_zinc.get('n_aizynth_solved', 0)} AiZynth / "
              f"{baseline_zinc.get('n_map_valid', 0)} map-valid / "
              f"{baseline_zinc.get('n_effective', 0)} 有效")
        print(f"  Strict:    {baseline_strict.get('n_aizynth_solved', 0)} AiZynth / "
              f"{baseline_strict.get('n_map_valid', 0)} map-valid / "
              f"{baseline_strict.get('n_effective', 0)} 有效")
        diff_eff = baseline_zinc.get('n_effective', 0) - baseline_strict.get('n_effective', 0)
        if diff_eff > 0:
            print(f"  → ZINC 比 Strict 多 {diff_eff} 条有效路线，"
                  f"说明有 {diff_eff} 条路线依赖复杂天然产物作为起始原料")
        elif diff_eff == 0:
            print(f"  → 两者有效路线数相同")

    # 2. Custom template 的影响 (ZINC stock)
    if baseline_zinc and flavonoid_zinc:
        print(f"\n[2] Custom template 影响 (ZINC stock):")
        print(f"  Baseline:  {baseline_zinc.get('n_effective', 0)} 有效")
        print(f"  +Custom:   {flavonoid_zinc.get('n_effective', 0)} 有效")
        diff_eff = flavonoid_zinc.get('n_effective', 0) - baseline_zinc.get('n_effective', 0)
        if diff_eff > 0:
            print(f"  → Custom template 增加了 {diff_eff} 条有效路线")
        elif diff_eff == 0:
            print(f"  → Custom template 未增加有效路线数")
        else:
            print(f"  → Custom template 减少了 {abs(diff_eff)} 条有效路线（可能引入噪声）")

    # 3. Custom template 的影响 (Strict stock)
    if baseline_strict and flavonoid_strict:
        print(f"\n[3] Custom template 影响 (Strict stock):")
        print(f"  Baseline+Strict: {baseline_strict.get('n_effective', 0)} 有效")
        print(f"  +Custom+Strict:  {flavonoid_strict.get('n_effective', 0)} 有效")
        diff_eff = flavonoid_strict.get('n_effective', 0) - baseline_strict.get('n_effective', 0)
        if diff_eff > 0:
            print(f"  → 在真实库存下，Custom template 增加了 {diff_eff} 条有效路线")
        elif diff_eff == 0:
            print(f"  → 在真实库存下，Custom template 未增加有效路线数")
        else:
            print(f"  → 在真实库存下，Custom template 减少了 {abs(diff_eff)} 条")

    # 4. Custom only 的贡献
    if custom_only:
        print(f"\n[4] 仅 Custom template (Strict stock):")
        print(f"  AiZynth solved: {custom_only.get('n_aizynth_solved', 0)}")
        print(f"  map-valid solved: {custom_only.get('n_map_valid', 0)}")
        print(f"  有效 solved: {custom_only.get('n_effective', 0)}")
        print(f"  自定义模板使用 (all routes): {custom_only.get('all_custom_family_counts', {})}")
        print(f"  到达 aglycone (all): {custom_only.get('reaches_aglycone_count_all', 0)}")
        print(f"  到达 chalcone (all): {custom_only.get('reaches_chalcone_count_all', 0)}")
        if custom_only.get('n_effective', 0) > 0:
            print(f"  → 黄酮模板自身可以独立找到有效解决路线")
        else:
            print(f"  → 黄酮模板自身无法独立找到有效解决路线，但在 unsolved branches 中确实进入了 MCTS")

    # 5. Sugar bridge closure 的影响
    if custom_only_virtual or flavonoid_virtual:
        print(f"\n[5] Sugar bridge closure 影响 (Virtual bridge):")
        if custom_only_virtual:
            print(
                f"  Custom only + virtual: effective={custom_only_virtual.get('n_effective', 0)}, "
                f"bridge-closed={custom_only_virtual.get('n_bridge_closed', 0)}, "
                f"non-virtual={custom_only_virtual.get('n_non_virtual_effective', 0)}"
            )
        if flavonoid_virtual:
            print(
                f"  USPTO+Custom + virtual: effective={flavonoid_virtual.get('n_effective', 0)}, "
                f"bridge-closed={flavonoid_virtual.get('n_bridge_closed', 0)}, "
                f"non-virtual={flavonoid_virtual.get('n_non_virtual_effective', 0)}"
            )
        print("  → 这些路线验证糖闭合层连通性，不计为 strict buyable solved")

    # 6. 零 retention 模板警告
    print(f"\n[6] 零 Map Retention 模板审计:")
    zero_ret_templates = defaultdict(list)
    for name, r in all_results.items():
        for a in r['analyses']:
            if a['map_retention'] == 0.0:
                for step in a['steps']:
                    tc = step.get('template_code', '')
                    cls = step.get('classification', '')
                    zero_ret_templates[(tc, cls)].append(name)

    if zero_ret_templates:
        print(f"  以下模板在所有实验中 map retention = 0（scaffold collapse）:")
        for (tc, cls), experiments in sorted(zero_ret_templates.items()):
            exps = ', '.join(set(experiments))
            print(f"  ⚠ template {tc} ({cls}): 出现在 {exps}")
    else:
        print(f"  无零 retention 模板")

    # 保存报告
    report_file = OUTPUT_DIR / "ablation_report.md"
    with open(report_file, 'w') as f:
        f.write("# 消融实验报告\n\n")
        f.write(f"生成日期: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")

        f.write("## 评价标准\n\n")
        f.write("- **AiZynth solved**: AiZynthFinder 原生 `is_solved` 字段\n")
        f.write("- **map-valid solved**: AiZynth solved 且所有反应步骤 atom-map retention ≥ 0.8\n")
        f.write("- **有效 solved**: map-valid solved（即真正通过守恒审计的路线）\n")
        f.write("- **bridge-closed solved**: 有效 solved 且至少一个叶子节点来自 `virtual_bridge`，只用于连通性验证，不等同于真实 strict solved\n")
        f.write("- **non-virtual effective solved**: 有效 solved 且不依赖 `virtual_bridge`\n")
        f.write("- **stock 层级 solved**: 路线中所有叶子节点均属于该 stock 层级\n")
        f.write("- **关键中间体到达**: 路线中包含到达 aglycone/chalcone/sugar donor 的步骤\n\n")

        f.write("## 汇总\n\n")
        f.write("| 实验 | 总路线 | AiZynth | map-valid | 有效 solved | bridge-closed | non-virtual | 最短步 |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|\n")
        for name, r in all_results.items():
            min_step = min(r['step_counts']) if r['step_counts'] else '-'
            f.write(f"| {name} | {r['n_total']} | {r['n_aizynth_solved']} | {r['n_map_valid']} | "
                    f"{r['n_effective']} | {r['n_bridge_closed']} | "
                    f"{r['n_non_virtual_effective']} | {min_step} |\n")

        f.write("\n## Stock 层级分析\n\n")
        f.write("| 实验 | strict_buyable | trusted_intermediate | virtual_bridge | unknown |\n")
        f.write("|---|---:|---:|---:|---:|\n")
        for name, r in all_results.items():
            layers = r['stock_layer_solved']
            f.write(f"| {name} | {layers.get('strict_buyable', 0)} | "
                    f"{layers.get('trusted_intermediate', 0)} | "
                    f"{layers.get('virtual_bridge', 0)} | "
                    f"{layers.get('unknown', 0)} |\n")

        f.write("\n## 自定义模板家族统计 (all routes)\n\n")
        f.write("| 实验 | 自定义模板家族 | 使用次数 |\n")
        f.write("|---|---|---:|\n")
        for name, r in all_results.items():
            if r['all_custom_family_counts']:
                for family, count in sorted(r['all_custom_family_counts'].items(), key=lambda x: -x[1]):
                    f.write(f"| {name} | {family} | {count} |\n")
            else:
                f.write(f"| {name} | - | 0 |\n")

        f.write("\n## 关键中间体到达统计\n\n")
        f.write("| 实验 | aglycone (solved/all) | chalcone (solved/all) | sugar donor (solved/all) |\n")
        f.write("|---|---:|---:|---:|\n")
        for name, r in all_results.items():
            f.write(f"| {name} | {r['reaches_aglycone_count']}/{r['reaches_aglycone_count_all']} | "
                    f"{r['reaches_chalcone_count']}/{r['reaches_chalcone_count_all']} | "
                    f"{r['reaches_sugar_donor_count']}/{r['reaches_sugar_donor_count_all']} |\n")

        f.write("\n## 结论\n\n")
        if baseline_zinc and baseline_strict:
            f.write(f"- Stock 影响 (Baseline): "
                    f"ZINC 有效={baseline_zinc.get('n_effective',0)}, "
                    f"Strict 有效={baseline_strict.get('n_effective',0)}\n")
        if baseline_strict and flavonoid_strict:
            f.write(f"- Custom template 影响 (Strict): "
                    f"Baseline 有效={baseline_strict.get('n_effective',0)}, "
                    f"+Custom 有效={flavonoid_strict.get('n_effective',0)}\n")
        if custom_only:
            f.write(f"- 仅 Custom template: "
                    f"AiZynth={custom_only.get('n_aizynth_solved',0)}, "
                    f"有效={custom_only.get('n_effective',0)}\n")
        if custom_only_virtual:
            f.write(f"- Sugar bridge (Custom only + virtual): "
                    f"有效={custom_only_virtual.get('n_effective',0)}, "
                    f"bridge-closed={custom_only_virtual.get('n_bridge_closed',0)}, "
                    f"non-virtual={custom_only_virtual.get('n_non_virtual_effective',0)}\n")
        if flavonoid_virtual:
            f.write(f"- Sugar bridge (USPTO+Custom + virtual): "
                    f"有效={flavonoid_virtual.get('n_effective',0)}, "
                    f"bridge-closed={flavonoid_virtual.get('n_bridge_closed',0)}, "
                    f"non-virtual={flavonoid_virtual.get('n_non_virtual_effective',0)}\n")

        # 零 retention 模板列表
        if zero_ret_templates:
            f.write("\n## 零 Map Retention 模板\n\n")
            f.write("| 模板代码 | 分类 | 出现实验 |\n")
            f.write("|---|---|---|\n")
            for (tc, cls), experiments in sorted(zero_ret_templates.items()):
                exps = ', '.join(set(experiments))
                f.write(f"| {tc} | {cls} | {exps} |\n")

    # 保存结构化 JSON 报告
    report_json = {
        'generated': __import__('datetime').datetime.now().isoformat(),
        'map_retention_threshold': MAP_RETENTION_THRESHOLD,
        'summary': {},
        'stock_layer_analysis': {},
        'custom_family_statistics': {},
        'intermediate_reaching_statistics': {},
        'zero_retention_templates': [
            {'template_code': tc, 'classification': cls, 'experiments': list(set(exps))}
            for (tc, cls), exps in sorted(zero_ret_templates.items())
        ],
    }
    for name, r in all_results.items():
        report_json['summary'][name] = {
            'n_total': r['n_total'],
            'n_aizynth_solved': r['n_aizynth_solved'],
            'n_map_valid': r['n_map_valid'],
            'n_effective': r['n_effective'],
            'n_bridge_closed': r['n_bridge_closed'],
            'n_non_virtual_effective': r['n_non_virtual_effective'],
            'min_steps': min(r['step_counts']) if r['step_counts'] else None,
            'custom_step_total_solved': r['custom_step_total'],
            'custom_step_total_all': r['custom_step_total_all'],
            'glyco_count_solved': r['glyco_count'],
            'glyco_count_all': r['glyco_count_all'],
            'chalcone_retro_aldol_count_solved': r['chalcone_retro_aldol_count'],
            'chalcone_retro_aldol_count_all': r['chalcone_retro_aldol_count_all'],
        }
        report_json['stock_layer_analysis'][name] = r['stock_layer_solved']
        report_json['custom_family_statistics'][name] = {
            'solved_only': r['custom_family_counts'],
            'all_routes': r['all_custom_family_counts'],
        }
        report_json['intermediate_reaching_statistics'][name] = {
            'reaches_aglycone_solved': r['reaches_aglycone_count'],
            'reaches_aglycone_all': r['reaches_aglycone_count_all'],
            'reaches_chalcone_solved': r['reaches_chalcone_count'],
            'reaches_chalcone_all': r['reaches_chalcone_count_all'],
            'reaches_sugar_donor_solved': r['reaches_sugar_donor_count'],
            'reaches_sugar_donor_all': r['reaches_sugar_donor_count_all'],
        }

    json_file = OUTPUT_DIR / "ablation_report.json"
    with open(json_file, 'w') as f:
        json.dump(report_json, f, indent=2, ensure_ascii=False)

    print(f"\n报告已保存: {report_file}")
    print(f"JSON 报告: {json_file}")


if __name__ == "__main__":
    compare_ablation()
