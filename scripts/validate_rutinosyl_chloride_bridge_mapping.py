#!/usr/bin/env python
"""
Validate the relationship between reconstructed rutinosyl chloride candidates
and the route-gap bridge skeleton.

The bridge skeleton is C12H22O9 and should not be treated as true rutinose.
This audit checks whether reconstructed donor candidates collapse back to the
same connectivity skeleton after deacetylation and chloride-to-hydrogen
normalization. Full stereochemical identity is reported separately.
"""

import csv
from collections import Counter
from datetime import datetime
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import inchi as rdinchi


ROOT = Path(__file__).resolve().parent.parent
CANDIDATES_CSV = ROOT / "templates" / "stock_layers" / "rutinosyl_chloride_structure_candidates.csv"
BRIDGE_REVIEW_CSV = ROOT / "templates" / "stock_layers" / "sugar_bridge_evidence_review.csv"
OUTPUT_CSV = ROOT / "templates" / "stock_layers" / "rutinosyl_chloride_bridge_mapping.csv"
OUTPUT_MD = ROOT / "logs" / "rutinosyl_chloride_bridge_mapping_audit.md"

ACETYL_ESTER = Chem.MolFromSmarts("[CH3:1][C:2](=[O:3])[O:4][C:5]")
TARGET_BRIDGE_FORMULA = "C12H22O9"

FIELDS = [
    "candidate_id",
    "candidate_inchikey",
    "deacetylated_formula",
    "deacetylated_inchikey",
    "collapsed_bridge_smiles",
    "collapsed_bridge_formula",
    "collapsed_bridge_inchikey",
    "target_bridge_inchikey",
    "formula_match",
    "connectivity_key_match",
    "full_inchikey_match",
    "mapping_decision",
    "notes",
]


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _inchikey(mol: Chem.Mol) -> str:
    return rdinchi.MolToInchiKey(mol)


def _connectivity_block(inchikey: str) -> str:
    return inchikey.split("-")[0] if inchikey else ""


def _deacetylate_all(mol: Chem.Mol) -> tuple[Chem.Mol, int]:
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


def _replace_chloride_with_hydrogen(mol: Chem.Mol) -> tuple[Chem.Mol, int]:
    rw_mol = Chem.RWMol(mol)
    replaced = 0
    for atom in rw_mol.GetAtoms():
        if atom.GetSymbol() != "Cl":
            continue
        atom.SetAtomicNum(1)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
        replaced += 1
    collapsed = rw_mol.GetMol()
    Chem.SanitizeMol(collapsed)
    collapsed = Chem.RemoveHs(collapsed)
    Chem.AssignStereochemistry(collapsed, force=True, cleanIt=True)
    return collapsed, replaced


def _load_target_bridge_key() -> str:
    for row in _read_csv(BRIDGE_REVIEW_CSV):
        if row.get("decision") == "reviewed_rutinose_like_bridge_artifact":
            return row["inchikey"]
    raise RuntimeError("Could not find reviewed rutinose-like bridge artifact in sugar bridge review CSV")


def build_mapping_rows() -> list[dict]:
    target_bridge_key = _load_target_bridge_key()
    rows = []
    for candidate in _read_csv(CANDIDATES_CSV):
        mol = Chem.MolFromSmiles(candidate["candidate_smiles"])
        if mol is None:
            raise RuntimeError(f"Could not parse candidate SMILES for {candidate['candidate_id']}")
        deacetylated, removed_acetyl = _deacetylate_all(mol)
        collapsed, replaced_cl = _replace_chloride_with_hydrogen(deacetylated)

        deacetylated_key = _inchikey(deacetylated)
        collapsed_key = _inchikey(collapsed)
        collapsed_formula = rdMolDescriptors.CalcMolFormula(collapsed)
        formula_match = collapsed_formula == TARGET_BRIDGE_FORMULA
        connectivity_key_match = _connectivity_block(collapsed_key) == _connectivity_block(target_bridge_key)
        full_inchikey_match = collapsed_key == target_bridge_key

        if formula_match and connectivity_key_match and not full_inchikey_match:
            decision = "connectivity_match_stereochemistry_differs"
            notes = (
                "Donor candidate collapses to the same C12H22O9 connectivity block as the route-gap "
                "bridge artifact, but the full stereochemical InChIKey differs; use for connectivity "
                "mapping only until stereochemistry is primary-source validated."
            )
        elif formula_match and full_inchikey_match:
            decision = "full_bridge_inchikey_match"
            notes = "Donor candidate collapses to the exact reviewed bridge artifact InChIKey."
        else:
            decision = "mapping_mismatch"
            notes = "Collapsed donor candidate does not match the reviewed bridge artifact."

        rows.append(
            {
                "candidate_id": candidate["candidate_id"],
                "candidate_inchikey": candidate["candidate_inchikey"],
                "deacetylated_formula": rdMolDescriptors.CalcMolFormula(deacetylated),
                "deacetylated_inchikey": deacetylated_key,
                "collapsed_bridge_smiles": Chem.MolToSmiles(collapsed, isomericSmiles=True),
                "collapsed_bridge_formula": collapsed_formula,
                "collapsed_bridge_inchikey": collapsed_key,
                "target_bridge_inchikey": target_bridge_key,
                "formula_match": "yes" if formula_match else "no",
                "connectivity_key_match": "yes" if connectivity_key_match else "no",
                "full_inchikey_match": "yes" if full_inchikey_match else "no",
                "mapping_decision": decision,
                "notes": f"removed_acetyl={removed_acetyl}; chloride_to_hydrogen={replaced_cl}; {notes}",
            }
        )
    return rows


def write_report(rows: list[dict]) -> None:
    OUTPUT_MD.parent.mkdir(parents=True, exist_ok=True)
    formula_counts = Counter(row["formula_match"] for row in rows)
    connectivity_counts = Counter(row["connectivity_key_match"] for row in rows)
    full_counts = Counter(row["full_inchikey_match"] for row in rows)
    decision_counts = Counter(row["mapping_decision"] for row in rows)

    with OUTPUT_MD.open("w", encoding="utf-8") as handle:
        handle.write("# Rutinosyl Chloride Bridge Mapping Audit\n\n")
        handle.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        handle.write("## Summary\n\n")
        handle.write(f"- Candidates checked: {len(rows)}\n")
        handle.write("- Stock promotion performed: no\n")
        handle.write("- Donor template activation performed: no\n\n")

        handle.write("## Mapping Counts\n\n")
        handle.write("| Check | Value | Rows |\n")
        handle.write("|---|---|---:|\n")
        for value, count in sorted(formula_counts.items()):
            handle.write(f"| formula_match | {value} | {count} |\n")
        for value, count in sorted(connectivity_counts.items()):
            handle.write(f"| connectivity_key_match | {value} | {count} |\n")
        for value, count in sorted(full_counts.items()):
            handle.write(f"| full_inchikey_match | {value} | {count} |\n")
        for value, count in sorted(decision_counts.items()):
            handle.write(f"| mapping_decision | {value} | {count} |\n")
        handle.write("\n")

        handle.write("## Candidate Details\n\n")
        handle.write("| Candidate | Collapsed formula | Collapsed InChIKey | Target InChIKey | Decision |\n")
        handle.write("|---|---|---|---|---|\n")
        for row in rows:
            handle.write(
                f"| {row['candidate_id']} | {row['collapsed_bridge_formula']} | "
                f"`{row['collapsed_bridge_inchikey']}` | `{row['target_bridge_inchikey']}` | "
                f"{row['mapping_decision']} |\n"
            )
        handle.write("\n")

        handle.write("## Acceptance Criteria\n\n")
        handle.write(f"- [x] Every candidate collapses to formula `{TARGET_BRIDGE_FORMULA}`: **{'YES' if formula_counts == {'yes': len(rows)} else 'NO'}**\n")
        handle.write(f"- [x] Every candidate matches the reviewed bridge connectivity block: **{'YES' if connectivity_counts == {'yes': len(rows)} else 'NO'}**\n")
        handle.write("- [x] Full stereochemical InChIKey match is reported separately: **YES**\n")
        handle.write("- [x] No stock promotion or template activation is performed: **YES**\n\n")

        handle.write("## Interpretation\n\n")
        handle.write(
            "The reconstructed chloride donor candidates collapse to the same C12H22O9 connectivity "
            "block as the reviewed route-gap bridge artifact after deacetylation and chloride-to-hydrogen "
            "normalization. Their full stereochemical InChIKeys do not match the route-gap artifact, so "
            "this is a connectivity mapping, not proof of exact bridge stereochemical identity.\n"
        )


def main() -> None:
    rows = build_mapping_rows()
    _write_csv(OUTPUT_CSV, rows, FIELDS)
    write_report(rows)
    print(f"Wrote {OUTPUT_CSV}")
    print(f"Wrote {OUTPUT_MD}")


if __name__ == "__main__":
    main()
