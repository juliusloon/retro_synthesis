#!/usr/bin/env python3
"""
构建 strict buyable stock 文件。

从 flavonoid_stock_inchikeys.txt 中过滤掉复杂黄酮糖苷和近邻天然产物，
只保留真正可购买的小分子前体（MW <= 350, 重原子数 <= 25）。

用法:
    python scripts/build_strict_stock.py [--mw-cutoff 350] [--heavy-atom-cutoff 25]
"""
import argparse
import csv
import sys
from datetime import datetime
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TEMPLATES_DIR = PROJECT_ROOT / "templates"
LOGS_DIR = PROJECT_ROOT / "logs"


def load_inchikeys(filepath: Path) -> list[str]:
    """读取 InChIKey 文件，每行一个。"""
    with open(filepath) as f:
        return [line.strip() for line in f if line.strip()]


def load_stock_candidates(csv_path: Path) -> dict[str, dict]:
    """从 stock_candidates.csv 加载 InChIKey -> {smiles, mw, heavy_atoms} 映射。"""
    result = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row["inchikey"]
            smi = row["smiles"]
            mw = float(row["mw"]) if row.get("mw") else None
            # 如果 CSV 中没有 mw，用 RDKit 计算
            if mw is None:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    mw = Descriptors.ExactMolWt(mol)
            # 计算重原子数
            mol = Chem.MolFromSmiles(smi)
            heavy_atoms = mol.GetNumHeavyAtoms() if mol else None
            result[key] = {"smiles": smi, "mw": mw, "heavy_atoms": heavy_atoms}
    return result


def build_strict_stock(
    source_keys: list[str],
    mol_info: dict[str, dict],
    mw_cutoff: float = 350.0,
    heavy_atom_cutoff: int = 25,
) -> tuple[list[str], list[dict], list[str]]:
    """过滤 stock InChIKeys，排除 MW 或重原子数超标的条目。

    返回 (kept, excluded, unverified)：
      - kept: 通过过滤的 InChIKey
      - excluded: 被排除的条目详情
      - unverified: 缺少分子信息的 InChIKey（不混入 strict stock）
    """
    kept = []
    excluded = []
    unverified = []

    for key in source_keys:
        if key not in mol_info:
            # 无法获取分子信息，不纳入 strict stock
            unverified.append(key)
            continue

        info = mol_info[key]
        mw = info["mw"]
        ha = info["heavy_atoms"]

        if mw is not None and mw > mw_cutoff:
            excluded.append({
                "inchikey": key,
                "smiles": info["smiles"],
                "mw": round(mw, 1),
                "heavy_atoms": ha,
                "reason": f"MW={round(mw,1)} > {mw_cutoff}",
            })
        elif ha is not None and ha > heavy_atom_cutoff:
            excluded.append({
                "inchikey": key,
                "smiles": info["smiles"],
                "mw": round(mw, 1) if mw else None,
                "heavy_atoms": ha,
                "reason": f"HA={ha} > {heavy_atom_cutoff}",
            })
        else:
            kept.append(key)

    return kept, excluded, unverified


def main():
    parser = argparse.ArgumentParser(description="构建 strict buyable stock 文件")
    parser.add_argument("--mw-cutoff", type=float, default=350.0, help="分子量上限（默认 350）")
    parser.add_argument("--heavy-atom-cutoff", type=int, default=25, help="重原子数上限（默认 25）")
    args = parser.parse_args()

    source_file = TEMPLATES_DIR / "flavonoid_stock_inchikeys.txt"
    candidates_csv = TEMPLATES_DIR / "literature_curated" / "flavonoid_literature_stock_candidates.csv"

    print(f"源 stock 文件: {source_file}")
    print(f"候选 CSV: {candidates_csv}")
    print(f"MW 截止: {args.mw_cutoff}")
    print(f"重原子数截止: {args.heavy_atom_cutoff}")
    print()

    # 加载数据
    source_keys = load_inchikeys(source_file)
    mol_info = load_stock_candidates(candidates_csv)
    print(f"源 stock 条目数: {len(source_keys)}")
    print(f"候选分子信息数: {len(mol_info)}")

    # 过滤
    strict_keys, excluded, unverified = build_strict_stock(
        source_keys, mol_info, args.mw_cutoff, args.heavy_atom_cutoff
    )

    # 输出 strict stock 文件
    output_file = TEMPLATES_DIR / "strict_buyable_stock_inchikeys.txt"
    with open(output_file, "w") as f:
        for key in strict_keys:
            f.write(key + "\n")

    # 输出 unverified stock 文件（单独存放，不混入 strict）
    unverified_file = TEMPLATES_DIR / "unverified_stock_inchikeys.txt"
    with open(unverified_file, "w") as f:
        for key in unverified:
            f.write(key + "\n")

    print(f"\n保留条目数: {len(strict_keys)}")
    print(f"排除条目数: {len(excluded)}")
    print(f"未验证条目数: {len(unverified)}")
    print(f"输出文件: {output_file}")
    print(f"未验证文件: {unverified_file}")

    # 生成审计报告
    report_file = LOGS_DIR / "26-06-29-strict_stock_audit.md"
    with open(report_file, "w") as f:
        f.write("# Strict Buyable Stock 审计报告\n\n")
        f.write(f"生成日期: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        f.write("## 参数\n\n")
        f.write(f"- 源文件: `{source_file.name}`\n")
        f.write(f"- MW 截止: {args.mw_cutoff}\n")
        f.write(f"- 重原子数截止: {args.heavy_atom_cutoff}\n\n")
        f.write("## 统计\n\n")
        f.write("| 指标 | 数值 |\n")
        f.write("|---|---:|\n")
        f.write(f"| 原始条目数 | {len(strict_keys) + len(excluded) + len(unverified)} |\n")
        f.write(f"| 保留条目数 | {len(strict_keys)} |\n")
        f.write(f"| 排除条目数 | {len(excluded)} |\n")
        f.write(f"| 未验证条目数 | {len(unverified)} |\n\n")

        f.write("## 处理策略\n\n")
        f.write("- **保留**: MW ≤ 350 且重原子数 ≤ 25 的可购买小分子\n")
        f.write("- **排除**: MW 或重原子数超标的复杂天然产物/糖苷\n")
        f.write("- **未验证**: 缺少分子信息的 InChIKey，**不纳入** strict stock（单独存放于 `unverified_stock_inchikeys.txt`）\n\n")

        if excluded:
            f.write("## 排除条目\n\n")
            f.write("| InChIKey | SMILES | MW | 重原子数 | 排除原因 |\n")
            f.write("|---|---|---:|---:|---|\n")
            for item in excluded:
                smi_short = item["smiles"][:60] + "..." if len(item["smiles"]) > 60 else item["smiles"]
                mw_str = str(item["mw"]) if item["mw"] else "N/A"
                ha_str = str(item["heavy_atoms"]) if item["heavy_atoms"] else "N/A"
                f.write(f"| `{item['inchikey']}` | `{smi_short}` | {mw_str} | {ha_str} | {item['reason']} |\n")
        else:
            f.write("## 无排除条目\n\n所有已验证 stock 条目均通过 strict 过滤。\n")

        if unverified:
            f.write(f"\n## 未验证条目 ({len(unverified)} 个)\n\n")
            f.write("以下 InChIKey 缺少分子信息（SMILES/MW），无法验证是否满足 strict 条件。\n")
            f.write("已单独存放于 `unverified_stock_inchikeys.txt`，不计入 strict stock。\n\n")
            for key in unverified[:20]:
                f.write(f"- `{key}`\n")
            if len(unverified) > 20:
                f.write(f"- ... 共 {len(unverified)} 个\n")

    print(f"审计报告: {report_file}")


if __name__ == "__main__":
    main()
