#!/usr/bin/env python3
"""
生成路线缺口诊断工作列表。

读取 unsolved routes，分析最常见未闭合的 leaf/intermediate，
输出工作列表以指导下一步的 stock 或 template 补充。

用法:
    python scripts/generate_route_gap_worklist.py
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, inchi

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "ablation"
STOCK_LAYERS_DIR = PROJECT_ROOT / "templates" / "stock_layers"


def load_all_stock_inchikeys() -> set:
    """加载所有 stock 层级的 InChIKey 集合。"""
    all_inchikeys = set()
    for layer_name in ["strict_buyable", "trusted_intermediate", "virtual_bridge"]:
        inchikey_file = STOCK_LAYERS_DIR / f"{layer_name}_stock_inchikeys.txt"
        if inchikey_file.exists():
            with open(inchikey_file) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        all_inchikeys.add(line)
    return all_inchikeys


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


def load_reaction_families() -> dict:
    """加载反应家族配置，返回分类到家族名称的映射。"""
    families_file = PROJECT_ROOT / "config" / "reaction_families.json"
    if not families_file.exists():
        return {}
    with open(families_file) as f:
        families = json.load(f)
    return {family["classification"]: family["name"] for family in families}


def smiles_to_inchikey(smiles: str) -> str:
    """从 SMILES 计算 InChIKey。"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        return inchi.MolToInchiKey(mol)
    except Exception:
        return ""


def canonicalize_smiles(smiles: str) -> str:
    """规范化 SMILES，用于路线节点匹配。"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return ""


def count_sugar_rings(mol) -> int:
    """统计非芳香、含氧、5/6元 C/O 糖环数量。"""
    count = 0
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) not in (5, 6):
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(atom.GetIsAromatic() for atom in atoms):
            continue
        symbols = [atom.GetSymbol() for atom in atoms]
        if "O" in symbols and all(symbol in {"C", "O"} for symbol in symbols):
            count += 1
    return count


def collect_unsolved_leaves(route_dict: dict) -> list:
    """收集路线中未解决的叶子节点（不在 stock 中的）。"""
    leaves = []

    def walk(node):
        if isinstance(node, dict):
            if node.get('type') == 'mol':
                children = node.get('children', [])
                if not children or len(children) == 0:
                    # 叶子节点
                    in_stock = node.get('in_stock', False)
                    if not in_stock:
                        smiles = node.get('smiles', '')
                        inchikey = node.get('inchikey', '')
                        # 如果没有 inchikey，从 SMILES 计算
                        if not inchikey and smiles:
                            inchikey = smiles_to_inchikey(smiles)
                        if smiles:
                            leaves.append({
                                'smiles': smiles,
                                'inchikey': inchikey,
                                'in_stock': in_stock
                            })
            # 所有类型节点都递归子节点
            for child in node.get('children', []):
                walk(child)

    walk(route_dict)
    return leaves


def get_upstream_family(route_dict: dict, target_smiles: str, target_inchikey: str = "") -> str:
    """获取导致生成目标 SMILES 的上游反应家族。"""
    reaction_families = load_reaction_families()
    target_canonical = canonicalize_smiles(target_smiles)
    target_inchikey = target_inchikey or smiles_to_inchikey(target_smiles)

    def child_matches_target(child):
        if child.get('type') != 'mol':
            return False
        child_smiles = child.get('smiles', '')
        child_inchikey = child.get('inchikey', '')
        if target_inchikey:
            child_inchikey = child_inchikey or smiles_to_inchikey(child_smiles)
            if child_inchikey == target_inchikey:
                return True
        return bool(target_canonical) and canonicalize_smiles(child_smiles) == target_canonical

    def walk(node, parent_family=None):
        if isinstance(node, dict):
            if node.get('type') == 'reaction':
                metadata = node.get('metadata', {})
                classification = metadata.get('classification', '')
                family_name = reaction_families.get(classification, classification)
                # 检查子节点是否包含目标 SMILES
                for child in node.get('children', []):
                    if child_matches_target(child):
                        return family_name
                    result = walk(child, family_name)
                    if result:
                        return result
            else:
                if node.get('type') == 'mol' and parent_family and child_matches_target(node):
                    return parent_family
                for child in node.get('children', []):
                    result = walk(child, parent_family)
                    if result:
                        return result
    return walk(route_dict)


def generate_gap_worklist():
    """生成路线缺口工作列表。"""
    # 加载所有 stock InChIKey
    all_stock_inchikeys = load_all_stock_inchikeys()
    stock_metadata = load_stock_metadata()

    # 收集所有 unsolved routes 中的叶子节点
    leaf_counts = defaultdict(lambda: {
        'count': 0,
        'first_seen_experiment': None,
        'upstream_family': None,
        'smiles': '',
        'inchikey': ''
    })

    result_files = sorted(OUTPUT_DIR.glob("*.json"))
    result_files = [f for f in result_files if f.name not in (
        "ablation_summary.json", "ablation_report.json", "route_gap_worklist.json")]

    if not result_files:
        print("错误: 未找到消融实验结果文件")
        print(f"请先运行: python scripts/run_ablation_experiments.py")
        sys.exit(1)

    print("分析 unsolved routes 中的叶子节点...")

    for filepath in result_files:
        exp_name = filepath.stem
        with open(filepath) as f:
            routes = json.load(f)

        for route in routes:
            # 只分析 unsolved routes
            is_solved = route.get('is_solved', False)
            if is_solved:
                continue

            # 收集未解决的叶子节点
            unsolved_leaves = collect_unsolved_leaves(route)

            for leaf in unsolved_leaves:
                smiles = leaf['smiles']
                inchikey = leaf['inchikey']

                if not smiles:
                    continue

                # 如果已经在 stock 中，跳过
                if inchikey and inchikey in all_stock_inchikeys:
                    continue

                # 更新计数（用 inchikey 或 smiles 作为 key）
                key = inchikey if inchikey else smiles
                leaf_counts[key]['count'] += 1
                leaf_counts[key]['smiles'] = smiles
                leaf_counts[key]['inchikey'] = inchikey

                if leaf_counts[key]['first_seen_experiment'] is None:
                    leaf_counts[key]['first_seen_experiment'] = exp_name

                # 获取上游反应家族
                if leaf_counts[key]['upstream_family'] is None:
                    upstream = get_upstream_family(route, smiles, inchikey)
                    if upstream:
                        leaf_counts[key]['upstream_family'] = upstream

    # 按 count 排序
    sorted_leaves = sorted(leaf_counts.values(), key=lambda x: -x['count'])

    # 生成建议操作
    worklist = []
    for leaf in sorted_leaves:
        smiles = leaf['smiles']
        inchikey = leaf['inchikey']
        count = leaf['count']

        # 分析分子特征，生成建议
        suggested_action = suggest_action(smiles, inchikey, count=count)

        worklist.append({
            'smiles': smiles,
            'inchikey': inchikey,
            'count': count,
            'first_seen_experiment': leaf['first_seen_experiment'],
            'upstream_family': leaf['upstream_family'] or 'unknown',
            'suggested_action': suggested_action,
            'notes': ''
        })

    # 输出 CSV
    output_file = OUTPUT_DIR / "route_gap_worklist.csv"
    with open(output_file, 'w', newline='') as f:
        fieldnames = ['smiles', 'inchikey', 'count', 'first_seen_experiment',
                      'upstream_family', 'suggested_action', 'notes']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(worklist)

    # 输出 JSON
    json_file = OUTPUT_DIR / "route_gap_worklist.json"
    with open(json_file, 'w') as f:
        json.dump(worklist, f, indent=2, ensure_ascii=False)

    print(f"\n生成工作列表完成")
    print(f"CSV 文件: {output_file}")
    print(f"JSON 文件: {json_file}")
    print(f"共 {len(worklist)} 个未解决的叶子节点")

    # 打印 Top 10
    print(f"\nTop 10 最常见未解决叶子节点:")
    print(f"{'排名':<4} {'计数':<6} {'SMILES':<50} {'建议操作'}")
    print("-" * 80)
    for i, item in enumerate(worklist[:10], 1):
        smiles_short = item['smiles'][:47] + "..." if len(item['smiles']) > 50 else item['smiles']
        print(f"{i:<4} {item['count']:<6} {smiles_short:<50} {item['suggested_action']}")


def suggest_action(smiles: str, inchikey: str, count: int = 0) -> str:
    """根据分子特征建议操作。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "reject_unreasonable"

    sugar_rings = count_sugar_rings(mol)
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O")
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")

    # 黄酮/芳香糖苷片段不能自动当 stock bridge，否则会掩盖核心骨架化学。
    if aromatic_atoms and sugar_rings:
        return "manual_review_aromatic_glycoside"

    # 非芳香糖片段由 sugar bridge layer 统一处理。只有重复出现的糖样叶子
    # 才自动进入 virtual_bridge；低频糖样叶子保留人工审查。
    if not aromatic_atoms and (sugar_rings or (oxygen_count >= 4 and carbon_count >= 5)):
        return "add_to_virtual_bridge" if count >= 2 else "manual_review_sugar_bridge"

    # 分子量检查
    mw = Descriptors.MolWt(mol)
    if mw > 500:
        return "check_commercial_availability"

    # 检查是否是常见 building block
    # 简单的启发式规则
    num_heavy_atoms = mol.GetNumHeavyAtoms()

    # 小分子（<=10 个重原子）可能是可购买的
    if num_heavy_atoms <= 10:
        return "add_to_trusted_stock"

    # 中等大小分子（11-20 个重原子）
    if num_heavy_atoms <= 20:
        # 检查是否是常见中间体
        if is_common_intermediate(mol):
            return "add_to_trusted_stock"
        else:
            return "check_commercial_availability"

    # 大分子（>20 个重原子）
    if num_heavy_atoms <= 30:
        return "add_to_virtual_bridge"

    # 非常大的分子
    return "reject_unreasonable"


def is_common_intermediate(mol) -> bool:
    """检查是否是常见中间体（简单的启发式规则）。"""
    # 检查是否含有常见官能团
    smarts_patterns = [
        "[OH]",  # 羟基
        "[CH3]",  # 甲基
        "[c]",  # 芳香碳
        "[C]=[O]",  # 羰基
        "[C][O][C]",  # 醚键
    ]

    match_count = 0
    for pattern in smarts_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            match_count += 1

    # 如果匹配多个常见官能团，可能是中间体
    return match_count >= 2


if __name__ == "__main__":
    generate_gap_worklist()
