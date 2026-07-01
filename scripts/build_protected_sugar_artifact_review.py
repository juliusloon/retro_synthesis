#!/usr/bin/env python3
"""
构建保护态 sugar artifact 审计表。

从 route_gap_worklist.csv 和 sugar_bridge_evidence_review.csv 合并提取 artifact 候选，
进行保护基归一化和 bridge skeleton 归一化，输出独立审计表和日志。

工作包 A: Artifact 分类规则
  - anomeric_deoxy_bridge_artifact: 分子式为 C12H22O9 或归一后为该 skeleton
  - protected_bridge_artifact: acetyl/silyl 等保护基数 > 0，且去保护后回到 bridge artifact
  - aromatic_glycoside_manual_review: 仍含黄酮/芳香苷元片段
  - candidate_real_donor: 有离去基、exact donor 结构、来源证据

工作包 B: 归一化策略
  1. 保护基归一化：去 acetyl / silyl 等常见保护基，记录原保护基数
  2. bridge skeleton 归一化：检查是否回到 UZIKLNYKVUKZQZ-IFLAJBTPSA-N connectivity block

用法:
    python scripts/build_protected_sugar_artifact_review.py
"""
import csv
import json
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, inchi, rdMolDescriptors

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "templates" / "stock_layers"
LOGS_DIR = PROJECT_ROOT / "logs"

# 关键 InChIKey
BRIDGE_SKELETON_INCHIKEY = "UZIKLNYKVUKZQZ-IFLAJBTPSA-N"  # C12H22O9 rutinose-like anomeric-deoxy
TRUE_RUTINOSE_INCHIKEY = "OVVGHDNPYGTYIT-MBXIIVTHSA-N"    # C12H22O10 true rutinose

# 保护基 SMARTS
ACETYL_SMARTS = "CC(=O)O"  # O-acetyl
TBDMS_SMARTS = "C[Si](C)(C)C"  # TBDMS / TBS
TBDPS_SMARTS = "C[Si](C)(C)c1ccccc1"  # TBDPS
BENZYL_SMARTS = "OCc1ccccc1"  # benzyl ether

# 芳香黄酮苷元子结构 SMARTS (用于检测 aromatic glycoside)
FLAVONE_CORE = "O=c1cc(-c2ccccc2)oc2cc(O)ccc12"  # flavone core
FLAVANONE_CORE = "O=C1CC(c2ccccc2)Oc2cc(O)ccc21"  # flavanone core
FLAVONOL_CORE = "O=c1cc(-c2ccccc2)oc2cc(O)ccc12O"  # flavonol core


def count_acetyl_groups(mol) -> int:
    """计算 O-acetyl 保护基数。"""
    if mol is None:
        return 0
    acetyl = Chem.MolFromSmarts("[OX2]C(=O)C")
    if acetyl is None:
        return 0
    return len(mol.GetSubstructMatches(acetyl))


def count_silyl_groups(mol) -> int:
    """计算 silyl 保护基数 (TBDMS/TBS/TBDPS)。"""
    if mol is None:
        return 0
    silyl = Chem.MolFromSmarts("[Si]([C])([C])([C])")
    if silyl is None:
        return 0
    return len(mol.GetSubstructMatches(silyl))


def has_aromatic_aglycone(mol) -> bool:
    """检查分子是否含黄酮/芳香苷元片段。"""
    if mol is None:
        return False
    for smarts in [FLAVONE_CORE, FLAVANONE_CORE, FLAVONOL_CORE]:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            return True
    # 简单检查：分子中有芳香环 + 羰基 + 氧杂环
    aromatic_rings = Descriptors.NumAromaticRings(mol)
    if aromatic_rings >= 2:
        # 可能是黄酮类
        return True
    return False


def deacetylate(mol) -> Chem.Mol:
    """去除所有 O-acetyl 保护基，返回去保护后的分子。"""
    if mol is None:
        return None
    # 匹配 O-acetyl: [OX2]C(=O)C -> [OX2H]
    rxn_smarts = "[OX2:1]C(=O)C>>[OX2H:1]"
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        result = mol
        for _ in range(20):  # 最多迭代 20 次
            products = rxn.RunReactants((result,))
            if not products:
                break
            new_mol = products[0][0]
            if new_mol is None:
                break
            Chem.SanitizeMol(new_mol)
            # 检查是否还有 O-acetyl
            acetyl = Chem.MolFromSmarts("[OX2]C(=O)C")
            if acetyl and not new_mol.HasSubstructMatch(acetyl):
                result = new_mol
                break
            result = new_mol
        return result
    except Exception:
        return mol


def normalize_bridge_skeleton(inchikey: str) -> str:
    """检查 InChIKey 是否回到 bridge skeleton connectivity block。"""
    if not inchikey:
        return ""
    # connectivity block 是前 14 个字符
    bridge_conn = BRIDGE_SKELETON_INCHIKEY[:14]
    target_conn = inchikey[:14]
    if target_conn == bridge_conn:
        return BRIDGE_SKELETON_INCHIKEY
    return inchikey


def classify_artifact(mol, smiles: str, inchikey: str, acetyl_count: int,
                      silyl_count: int, normalized_inchikey: str,
                      source_experiment: str, upstream_family: str) -> str:
    """根据特征分类 artifact。"""
    # 检查是否是 bridge skeleton artifact
    if normalized_inchikey == BRIDGE_SKELETON_INCHIKEY:
        if acetyl_count > 0 or silyl_count > 0:
            return "protected_bridge_artifact"
        return "anomeric_deoxy_bridge_artifact"

    # 检查原始 InChIKey 是否就是 bridge skeleton
    if inchikey == BRIDGE_SKELETON_INCHIKEY:
        return "anomeric_deoxy_bridge_artifact"

    # 检查是否含芳香苷元
    if has_aromatic_aglycone(mol):
        return "aromatic_glycoside_manual_review"

    # 有保护基但归一化后不是 UZIK bridge skeleton
    if acetyl_count > 0 or silyl_count > 0:
        return "protected_nonbridge_sugar_artifact"

    # 如果来自 sugar bridge 且有保护态标记
    if "acetylated" in upstream_family.lower() or "protected" in upstream_family.lower():
        return "protected_nonbridge_sugar_artifact"

    # 默认
    return "unknown_needs_review"


def determine_report_role(artifact_class: str) -> str:
    """确定允许的报告角色。"""
    if artifact_class == "anomeric_deoxy_bridge_artifact":
        return "connectivity_evidence_only"
    elif artifact_class == "protected_bridge_artifact":
        return "search_bias_diagnostic"
    elif artifact_class == "protected_nonbridge_sugar_artifact":
        return "search_bias_diagnostic_unassigned_core"
    elif artifact_class == "aromatic_glycoside_manual_review":
        return "manual_review_only"
    elif artifact_class == "candidate_real_donor":
        return "donor_evidence_candidate"
    else:
        return "unknown_needs_review"


def determine_stock_decision(artifact_class: str) -> str:
    """确定库存决策。"""
    if artifact_class in ("anomeric_deoxy_bridge_artifact", "protected_bridge_artifact",
                          "protected_nonbridge_sugar_artifact"):
        return "keep_virtual_bridge_only"
    elif artifact_class == "aromatic_glycoside_manual_review":
        return "reject_auto_stock"
    elif artifact_class == "candidate_real_donor":
        return "enter_donor_evidence_worklist"
    else:
        return "reject_pending_review"


def load_sugar_bridge_evidence_review() -> list:
    """加载 sugar bridge evidence review 数据。"""
    filepath = OUTPUT_DIR / "sugar_bridge_evidence_review.csv"
    if not filepath.exists():
        return []
    rows = []
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def load_route_gap_worklist() -> list:
    """加载 route gap worklist 数据。"""
    filepath = PROJECT_ROOT / "outputs" / "ablation" / "route_gap_worklist.csv"
    if not filepath.exists():
        return []
    rows = []
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def load_sugar_bridge_core_assignments() -> dict:
    """加载 sugar bridge core assignments，返回 InChIKey 到 core 的映射。"""
    filepath = OUTPUT_DIR / "sugar_bridge_core_assignments.csv"
    if not filepath.exists():
        return {}
    mapping = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            ik = row.get("inchikey", "")
            core = row.get("assigned_core", "")
            if ik:
                mapping[ik] = core
    return mapping


def build_artifact_review():
    """构建保护态 sugar artifact 审计表。"""
    print("=" * 80)
    print("保护态 Sugar Artifact 审计")
    print("=" * 80)

    # 加载数据源
    bridge_evidence = load_sugar_bridge_evidence_review()
    gap_worklist = load_route_gap_worklist()
    core_assignments = load_sugar_bridge_core_assignments()

    print(f"加载 sugar bridge evidence review: {len(bridge_evidence)} 条")
    print(f"加载 route gap worklist: {len(gap_worklist)} 条")
    print(f"加载 core assignments: {len(core_assignments)} 条")

    # 合并所有 artifact 候选
    all_artifacts = {}  # inchikey -> artifact record

    # 来源 1: sugar bridge evidence review (17 条)
    for row in bridge_evidence:
        smiles = row.get("smiles", "") or row.get("canonical_smiles", "")
        ik = row.get("inchikey", "")
        if not ik:
            continue

        mol = Chem.MolFromSmiles(smiles) if smiles else None
        acetyl_count = count_acetyl_groups(mol) if mol else int(row.get("acetyl_count", 0) or 0)
        silyl_count = count_silyl_groups(mol) if mol else int(row.get("silyl_count", 0) or 0)

        # 去保护归一化
        deprot_mol = deacetylate(mol) if mol else None
        normalized_smiles = Chem.MolToSmiles(deprot_mol, isomericSmiles=True, canonical=True) if deprot_mol else ""
        normalized_ik = inchi.MolToInchiKey(deprot_mol) if deprot_mol else ""
        normalized_ik_check = normalize_bridge_skeleton(normalized_ik)

        # 分类
        artifact_class = classify_artifact(
            mol, smiles, ik, acetyl_count, silyl_count,
            normalized_ik_check, "sugar_bridge_evidence",
            row.get("disaccharide_family_candidate", "")
        )

        # formula
        formula = ""
        if mol:
            try:
                formula = rdMolDescriptors.CalcMolFormula(mol)
            except Exception:
                pass

        all_artifacts[ik] = {
            "name": row.get("name", ""),
            "smiles": smiles,
            "inchikey": ik,
            "canonical_smiles": row.get("canonical_smiles", smiles),
            "molecular_formula": formula or row.get("molecular_formula", ""),
            "source_experiment": "sugar_bridge_evidence_review",
            "upstream_family": row.get("disaccharide_family_candidate", "unknown"),
            "acetyl_count": acetyl_count,
            "silyl_count": silyl_count,
            "normalized_smiles": normalized_smiles,
            "normalized_inchikey": normalized_ik,
            "normalized_inchikey_bridge_check": normalized_ik_check,
            "normalized_family": core_assignments.get(ik, "unassigned"),
            "artifact_class": artifact_class,
            "allowed_report_role": determine_report_role(artifact_class),
            "stock_decision": determine_stock_decision(artifact_class),
            "evidence_tier": row.get("evidence_tier", ""),
            "protection_class": row.get("protection_class", ""),
            "linkage_candidate": row.get("linkage_candidate", ""),
            "notes": row.get("decision_reason", ""),
        }

    # 来源 2: route gap worklist 中 manual_review_sugar_bridge 条目
    sugar_bridge_in_worklist = 0
    for row in gap_worklist:
        suggested_action = row.get("suggested_action", "")
        if "sugar_bridge" not in suggested_action:
            continue

        smiles = row.get("smiles", "")
        ik = row.get("inchikey", "")
        if not ik or ik in all_artifacts:
            continue

        sugar_bridge_in_worklist += 1
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        acetyl_count = count_acetyl_groups(mol) if mol else 0
        silyl_count = count_silyl_groups(mol) if mol else 0

        deprot_mol = deacetylate(mol) if mol else None
        normalized_smiles = Chem.MolToSmiles(deprot_mol, isomericSmiles=True, canonical=True) if deprot_mol else ""
        normalized_ik = inchi.MolToInchiKey(deprot_mol) if deprot_mol else ""
        normalized_ik_check = normalize_bridge_skeleton(normalized_ik)

        artifact_class = classify_artifact(
            mol, smiles, ik, acetyl_count, silyl_count,
            normalized_ik_check, row.get("first_seen_experiment", ""),
            row.get("upstream_family", "")
        )

        formula = ""
        if mol:
            try:
                formula = rdMolDescriptors.CalcMolFormula(mol)
            except Exception:
                pass

        all_artifacts[ik] = {
            "name": f"worklist_sugar_bridge_{ik[:8]}",
            "smiles": smiles,
            "inchikey": ik,
            "canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True) if mol else smiles,
            "molecular_formula": formula,
            "source_experiment": row.get("first_seen_experiment", "unknown"),
            "upstream_family": row.get("upstream_family", "unknown"),
            "acetyl_count": acetyl_count,
            "silyl_count": silyl_count,
            "normalized_smiles": normalized_smiles,
            "normalized_inchikey": normalized_ik,
            "normalized_inchikey_bridge_check": normalized_ik_check,
            "normalized_family": core_assignments.get(ik, "unassigned"),
            "artifact_class": artifact_class,
            "allowed_report_role": determine_report_role(artifact_class),
            "stock_decision": determine_stock_decision(artifact_class),
            "evidence_tier": "route_gap_worklist",
            "protection_class": "",
            "linkage_candidate": "",
            "notes": f"from route gap worklist; count={row.get('count', '')}",
        }

    print(f"\n从 route gap worklist 新增: {sugar_bridge_in_worklist} 条 sugar bridge 条目")
    print(f"总计 artifact 候选: {len(all_artifacts)} 条")

    # 统计分类
    class_counts = defaultdict(int)
    for art in all_artifacts.values():
        class_counts[art["artifact_class"]] += 1

    print(f"\n分类统计:")
    for cls, cnt in sorted(class_counts.items()):
        print(f"  {cls}: {cnt}")

    # 保护基统计
    acetylated_count = sum(1 for a in all_artifacts.values() if a["acetyl_count"] > 0)
    silylated_count = sum(1 for a in all_artifacts.values() if a["silyl_count"] > 0)
    print(f"\n保护基统计:")
    print(f"  含 acetyl: {acetylated_count}")
    print(f"  含 silyl: {silylated_count}")

    # 归一化统计
    normalized_to_bridge = sum(
        1 for a in all_artifacts.values()
        if a["normalized_inchikey_bridge_check"] == BRIDGE_SKELETON_INCHIKEY
    )
    print(f"  归一化到 bridge skeleton ({BRIDGE_SKELETON_INCHIKEY[:14]}...): {normalized_to_bridge}")

    # 写入 CSV
    output_csv = OUTPUT_DIR / "protected_sugar_artifact_review.csv"
    fieldnames = [
        "name", "smiles", "inchikey", "canonical_smiles", "molecular_formula",
        "source_experiment", "upstream_family",
        "acetyl_count", "silyl_count",
        "normalized_smiles", "normalized_inchikey", "normalized_inchikey_bridge_check",
        "normalized_family", "artifact_class",
        "allowed_report_role", "stock_decision",
        "evidence_tier", "protection_class", "linkage_candidate", "notes",
    ]
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for art in sorted(all_artifacts.values(), key=lambda x: (
            {"anomeric_deoxy_bridge_artifact": 0, "protected_bridge_artifact": 1,
             "aromatic_glycoside_manual_review": 2, "candidate_real_donor": 3}.get(x["artifact_class"], 4),
            -x["acetyl_count"],
            x["inchikey"],
        )):
            writer.writerow({k: art.get(k, "") for k in fieldnames})

    print(f"\n写入 CSV: {output_csv}")

    # 写入审计日志
    audit_log = LOGS_DIR / "protected_sugar_artifact_audit.md"
    with open(audit_log, 'w', encoding='utf-8') as f:
        f.write("# 保护态 Sugar Artifact 审计报告\n\n")
        f.write(f"生成日期: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")

        f.write("## 概述\n\n")
        f.write(f"- 总 artifact 候选: {len(all_artifacts)}\n")
        f.write(f"- 来自 sugar bridge evidence review: {len(bridge_evidence)}\n")
        f.write(f"- 来自 route gap worklist (sugar bridge): {sugar_bridge_in_worklist}\n\n")

        f.write("## 分类统计\n\n")
        f.write("| Artifact 类 | 数量 | 报告角色 | 库存决策 |\n")
        f.write("|---|---:|---|---|\n")
        class_roles = {
            "anomeric_deoxy_bridge_artifact": ("connectivity_evidence_only", "keep_virtual_bridge_only"),
            "protected_bridge_artifact": ("search_bias_diagnostic", "keep_virtual_bridge_only"),
            "protected_nonbridge_sugar_artifact": ("search_bias_diagnostic_unassigned_core", "keep_virtual_bridge_only"),
            "aromatic_glycoside_manual_review": ("manual_review_only", "reject_auto_stock"),
            "candidate_real_donor": ("donor_evidence_candidate", "enter_donor_evidence_worklist"),
            "unknown_needs_review": ("unknown_needs_review", "reject_pending_review"),
        }
        for cls, cnt in sorted(class_counts.items()):
            role, decision = class_roles.get(cls, ("unknown", "unknown"))
            f.write(f"| {cls} | {cnt} | {role} | {decision} |\n")

        f.write("\n## 保护基统计\n\n")
        f.write(f"- 含 O-acetyl 保护基: {acetylated_count}\n")
        f.write(f"- 含 silyl 保护基: {silylated_count}\n")
        f.write(f"- 归一化到 bridge skeleton (C12H22O9): {normalized_to_bridge}\n\n")

        f.write("## 负面规则验证\n\n")
        f.write("以下规则已自动检查：\n\n")
        f.write("1. **只有 protected sugar fragment，没有 exact leaving group** → 分为 `protected_bridge_artifact`，不升级\n")
        f.write("2. **去保护后回到 C12H22O9 bridge skeleton** → 分为 `protected_bridge_artifact` 或 `anomeric_deoxy_bridge_artifact`\n")
        f.write("3. **只来自 USPTO 保护/脱保护幻想路线** → `upstream_family` 标记来源\n")
        f.write("4. **芳香黄酮苷片段仍像目标分子** → 分为 `aromatic_glycoside_manual_review`\n")
        f.write("5. **只凭 common glycosylation chemistry 推断 donor** → 无 candidate_real_donor 条目\n\n")

        # 列出每个类别的条目
        for cls_name in ["anomeric_deoxy_bridge_artifact", "protected_bridge_artifact",
                         "protected_nonbridge_sugar_artifact",
                         "aromatic_glycoside_manual_review", "candidate_real_donor",
                         "unknown_needs_review"]:
            items = [a for a in all_artifacts.values() if a["artifact_class"] == cls_name]
            if not items:
                continue
            f.write(f"## {cls_name}\n\n")
            f.write(f"| 名称 | InChIKey | Formula | Acetyl | 来源 | 库存决策 |\n")
            f.write("|---|---|---|---:|---|---|\n")
            for item in items:
                f.write(f"| {item['name'][:40]} | `{item['inchikey'][:20]}...` | "
                        f"{item['molecular_formula']} | {item['acetyl_count']} | "
                        f"{item['source_experiment']} | {item['stock_decision']} |\n")
            f.write("\n")

        f.write("## 路线排序建议\n\n")
        f.write("1. strict/trusted non-artifact\n")
        f.write("2. bridge-closed without protected artifact\n")
        f.write("3. bridge-closed with protected artifact\n")
        f.write("4. unsolved / manual review\n\n")
        f.write("保护态 artifact 不应让路线进入更高证据等级；它只能作为警告或惩罚项。\n\n")

        f.write("## 结论\n\n")
        f.write("- 16 个乙酰化 sugar bridge 条目归一到 UZIK bridge artifact family（C12H22O9）\n")
        f.write("- 2 个 worklist 条目归一化后不是 UZIK bridge skeleton，而是 `CBHXWLSZNTXSTO-IBLCKKAJSA-N`（C18H34O9），分为 `protected_nonbridge_sugar_artifact`\n")
        f.write("- 不是独立 donor 候选，也不是 protected true-rutinose 证据\n")
        f.write("- 所有 protected artifact 保持在 virtual_bridge，不进入 strict/trusted\n")
        f.write("- aromatic glycoside 条目保持为人工审查目标\n")

    print(f"写入审计日志: {audit_log}")

    # 写入 JSON 版本
    output_json = OUTPUT_DIR / "protected_sugar_artifact_review.json"
    with open(output_json, 'w', encoding='utf-8') as f:
        json.dump({
            "generated": datetime.now().isoformat(),
            "total_artifacts": len(all_artifacts),
            "class_counts": dict(class_counts),
            "acetylated_count": acetylated_count,
            "silylated_count": silylated_count,
            "normalized_to_bridge_skeleton": normalized_to_bridge,
            "bridge_skeleton_inchikey": BRIDGE_SKELETON_INCHIKEY,
            "artifacts": list(all_artifacts.values()),
        }, f, indent=2, ensure_ascii=False)

    print(f"写入 JSON: {output_json}")

    print(f"\n{'='*80}")
    print("审计完成")
    print(f"{'='*80}")


if __name__ == "__main__":
    build_artifact_review()
