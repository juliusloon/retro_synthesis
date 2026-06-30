#!/usr/bin/env python3
"""Group protected sugar bridge entries by deprotected sugar core.

This script is intentionally conservative: it does not promote any stock layer.
It removes simple O-acetyl groups from route-gap sugar bridge entries, compares
the resulting free-sugar core to manually reviewed core rows, and emits a
review worklist for only the entries that fail assignment.
"""

import csv
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import inchi

PROJECT_ROOT = Path(__file__).resolve().parent.parent
STOCK_CSV = PROJECT_ROOT / "templates" / "stock_layers" / "sugar_bridge_stock.csv"
EVIDENCE_CSV = PROJECT_ROOT / "templates" / "stock_layers" / "sugar_bridge_evidence_review.csv"
OUTPUT_CSV = PROJECT_ROOT / "templates" / "stock_layers" / "sugar_bridge_core_assignments.csv"
OUTPUT_MD = PROJECT_ROOT / "logs" / "sugar_bridge_core_assignment.md"

FIELDS = [
    "name",
    "inchikey",
    "protection_class",
    "acetyl_count",
    "deacetylated_groups_removed",
    "normalized_core_smiles",
    "normalized_core_inchikey",
    "normalized_core_assignment",
    "assignment_confidence",
    "human_review_needed",
    "decision",
    "notes",
]

ACETYL_ESTER = Chem.MolFromSmarts("[CH3:1][C:2](=[O:3])[O:4][C:5]")


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fileobj:
        return list(csv.DictReader(fileobj))


def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fileobj:
        writer = csv.DictWriter(fileobj, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _mol_from_smiles(smiles: str) -> Chem.Mol | None:
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def _inchikey(mol: Chem.Mol | None) -> str:
    if mol is None:
        return ""
    try:
        return inchi.MolToInchiKey(mol)
    except Exception:
        return ""


def _smiles(mol: Chem.Mol | None) -> str:
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def deacetylate_all(mol: Chem.Mol | None) -> tuple[Chem.Mol | None, int]:
    """Remove all simple O-acetyl groups while preserving sugar O atoms."""
    if mol is None:
        return None, 0
    current = Chem.Mol(mol)
    removed = 0
    while True:
        matches = current.GetSubstructMatches(ACETYL_ESTER)
        if not matches:
            return current, removed
        methyl, carbonyl, oxo, _ester_o, _sugar_c = matches[0]
        rw_mol = Chem.RWMol(current)
        for atom_idx in sorted([methyl, carbonyl, oxo], reverse=True):
            rw_mol.RemoveAtom(atom_idx)
        current = rw_mol.GetMol()
        Chem.SanitizeMol(current)
        removed += 1
        if removed > 20:
            raise RuntimeError("Too many acetyl removals; possible malformed molecule")


def _load_reviewed_cores(review_rows: list[dict]) -> dict[str, str]:
    """Return normalized core InChIKey to assignment label."""
    cores = {}
    for row in review_rows:
        decision = row.get("decision", "")
        linkage = row.get("linkage_candidate", "")
        inchikey_value = row.get("inchikey", "")
        if decision == "normalize_to_known_sugar_core" and "rutinose" in linkage:
            cores[inchikey_value] = "rutinose_core"
        elif decision == "reviewed_rutinose_like_bridge_artifact" and "rutinose" in linkage:
            cores[inchikey_value] = "rutinose_like_bridge_core_missing_anomeric_oxygen"
    return cores


def build_core_assignments() -> None:
    stock_rows = _read_csv(STOCK_CSV)
    review_rows = _read_csv(EVIDENCE_CSV)
    review_by_key = {row["inchikey"]: row for row in review_rows}
    reviewed_cores = _load_reviewed_cores(review_rows)

    assignment_rows = []
    for stock_row in stock_rows:
        inchikey_value = stock_row.get("inchikey", "")
        review_row = review_by_key.get(inchikey_value, {})
        mol = _mol_from_smiles(stock_row.get("smiles", ""))
        core_mol, removed = deacetylate_all(mol)
        core_key = _inchikey(core_mol)
        core_smiles = _smiles(core_mol)
        assignment = reviewed_cores.get(core_key, "")

        if assignment:
            confidence = "high"
            human_review_needed = "no"
            decision = "derived_from_reviewed_core"
            notes = (
                f"deacetylation maps this bridge to {assignment}; "
                "not independent evidence and not promoted from virtual_bridge"
            )
        else:
            confidence = "low"
            human_review_needed = "yes"
            decision = "needs_human_review"
            notes = "deacetylated core does not match a manually reviewed core"

        assignment_rows.append(
            {
                "name": stock_row.get("name", ""),
                "inchikey": inchikey_value,
                "protection_class": review_row.get("protection_class", ""),
                "acetyl_count": review_row.get("acetyl_count", ""),
                "deacetylated_groups_removed": removed,
                "normalized_core_smiles": core_smiles,
                "normalized_core_inchikey": core_key,
                "normalized_core_assignment": assignment,
                "assignment_confidence": confidence,
                "human_review_needed": human_review_needed,
                "decision": decision,
                "notes": notes,
            }
        )

    _write_csv(OUTPUT_CSV, assignment_rows, FIELDS)
    _write_report(assignment_rows, reviewed_cores)

    print(f"Wrote sugar bridge core assignments: {OUTPUT_CSV}")
    print(f"Wrote sugar bridge core assignment report: {OUTPUT_MD}")


def _write_report(rows: list[dict], reviewed_cores: dict[str, str]) -> None:
    OUTPUT_MD.parent.mkdir(parents=True, exist_ok=True)
    assignment_counts = Counter(row["normalized_core_assignment"] or "unassigned" for row in rows)
    review_counts = Counter(row["human_review_needed"] for row in rows)
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["normalized_core_assignment"] or "unassigned"].append(row)

    with OUTPUT_MD.open("w", encoding="utf-8") as fileobj:
        fileobj.write("# Sugar Bridge Core Assignment\n\n")
        fileobj.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        fileobj.write("## Summary\n\n")
        fileobj.write(f"- Sugar bridge entries assigned: {len(rows)}\n")
        fileobj.write(f"- Reviewed cores: {len(reviewed_cores)}\n")
        fileobj.write(f"- Human-review-needed entries: {review_counts.get('yes', 0)}\n")
        fileobj.write("- Stock promotion performed: no\n\n")

        fileobj.write("## Assignment Counts\n\n")
        fileobj.write("| Assignment | Rows |\n")
        fileobj.write("|---|---:|\n")
        for assignment, count in assignment_counts.most_common():
            fileobj.write(f"| {assignment} | {count} |\n")
        fileobj.write("\n")

        fileobj.write("## Human Review Counts\n\n")
        fileobj.write("| Human review needed | Rows |\n")
        fileobj.write("|---|---:|\n")
        for value, count in review_counts.most_common():
            fileobj.write(f"| {value} | {count} |\n")
        fileobj.write("\n")

        fileobj.write("## Assignment Details\n\n")
        fileobj.write("| Name | InChIKey | Removed acetyl groups | Core InChIKey | Assignment | Review needed |\n")
        fileobj.write("|---|---|---:|---|---|---|\n")
        for row in rows:
            fileobj.write(
                f"| {row['name']} | `{row['inchikey']}` | "
                f"{row['deacetylated_groups_removed']} | "
                f"`{row['normalized_core_inchikey']}` | "
                f"{row['normalized_core_assignment'] or 'unassigned'} | "
                f"{row['human_review_needed']} |\n"
            )
        fileobj.write("\n")

        fileobj.write("## Interpretation\n\n")
        if review_counts.get("yes", 0) == 0:
            fileobj.write(
                "All current sugar bridge entries deacetylate to a manually reviewed "
                "rutinose-like bridge skeleton. The reviewed core has formula "
                "C12H22O9, not true rutinose C12H22O10, so treat protected entries "
                "as one anomeric-deoxy bridge artifact family rather than independent "
                "trusted/strict evidence or free-rutinose donor evidence.\n"
            )
        else:
            fileobj.write(
                "Some entries do not map to a reviewed core and should be inspected "
                "before changing evidence tiers or donor design assumptions.\n"
            )


def main() -> None:
    build_core_assignments()


if __name__ == "__main__":
    main()
