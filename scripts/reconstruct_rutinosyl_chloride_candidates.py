#!/usr/bin/env python
"""
Reconstruct machine-checkable hexa-O-acetyl rutinosyl chloride candidates.

The literature hit names hexa-O-acetyl-beta-rutinosyl chloride, but the primary
paper structure is not machine-readable in the local pipeline. This script
builds conservative candidate anomers from the curated true-rutinose cyclic
seed, then records them as reconstruction candidates only.

It does not promote stock and does not enable donor templates.
"""

import csv
from collections import Counter
from datetime import datetime
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import inchi as rdinchi
from rdkit.Chem.rdchem import ChiralType


ROOT = Path(__file__).resolve().parent.parent
METADATA_CSV = ROOT / "templates" / "stock_layers" / "stock_layers_metadata.csv"
OUTPUT_CSV = ROOT / "templates" / "stock_layers" / "rutinosyl_chloride_structure_candidates.csv"
OUTPUT_MD = ROOT / "logs" / "rutinosyl_chloride_reconstruction_audit.md"

TARGET_NAME = "hexa-O-acetyl-beta-rutinosyl chloride"
TARGET_DOI = "https://doi.org/10.1016/S0008-6215(00)84374-8"
EXPECTED_FORMULA = "C24H33ClO15"

FIELDS = [
    "candidate_id",
    "literature_name",
    "reconstruction_label",
    "source_seed_name",
    "source_seed_inchikey",
    "source_seed_smiles",
    "anomeric_chiral_tag",
    "anomeric_assignment",
    "candidate_smiles",
    "candidate_inchikey",
    "candidate_formula",
    "exact_mass",
    "acetyl_count",
    "residual_oh_count",
    "chlorine_count",
    "formula_check",
    "evidence_status",
    "stock_decision",
    "template_decision",
    "validation_notes",
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


def _load_rutinose_seed() -> dict:
    for row in _read_csv(METADATA_CSV):
        if row.get("name") == "rutinose":
            return row
    raise RuntimeError("Could not find rutinose seed in stock_layers_metadata.csv")


def _terminal_oh_oxygen_indices(mol: Chem.Mol) -> list[int]:
    indices = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        if atom.GetDegree() == 1 and atom.GetTotalNumHs() == 1:
            indices.append(atom.GetIdx())
    return indices


def _find_reducing_anomeric_oh(mol: Chem.Mol) -> tuple[int, int]:
    candidates = []
    for oxygen_idx in _terminal_oh_oxygen_indices(mol):
        oxygen = mol.GetAtomWithIdx(oxygen_idx)
        carbon = oxygen.GetNeighbors()[0]
        ring_oxygen_neighbors = [
            atom for atom in carbon.GetNeighbors()
            if atom.GetIdx() != oxygen_idx and atom.GetSymbol() == "O" and atom.IsInRing()
        ]
        if carbon.IsInRing() and ring_oxygen_neighbors:
            candidates.append((oxygen_idx, carbon.GetIdx()))
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one reducing-end anomeric OH, found {len(candidates)}")
    return candidates[0]


def _acetylate_oxygen(rw_mol: Chem.RWMol, oxygen_idx: int) -> None:
    carbonyl = Chem.Atom("C")
    carbonyl.SetNoImplicit(True)
    carbonyl_idx = rw_mol.AddAtom(carbonyl)

    oxo = Chem.Atom("O")
    oxo.SetNoImplicit(True)
    oxo_idx = rw_mol.AddAtom(oxo)

    methyl = Chem.Atom("C")
    methyl_idx = rw_mol.AddAtom(methyl)

    rw_mol.AddBond(oxygen_idx, carbonyl_idx, Chem.BondType.SINGLE)
    rw_mol.AddBond(carbonyl_idx, oxo_idx, Chem.BondType.DOUBLE)
    rw_mol.AddBond(carbonyl_idx, methyl_idx, Chem.BondType.SINGLE)


def _make_candidate(seed_mol: Chem.Mol, anomeric_oh_idx: int, anomeric_c_idx: int, tag: ChiralType) -> Chem.Mol:
    rw_mol = Chem.RWMol(seed_mol)

    leaving_group = rw_mol.GetAtomWithIdx(anomeric_oh_idx)
    leaving_group.SetAtomicNum(17)
    leaving_group.SetFormalCharge(0)
    leaving_group.SetNoImplicit(True)
    leaving_group.SetNumExplicitHs(0)

    anomeric_carbon = rw_mol.GetAtomWithIdx(anomeric_c_idx)
    anomeric_carbon.SetChiralTag(tag)

    for oxygen_idx in _terminal_oh_oxygen_indices(rw_mol):
        _acetylate_oxygen(rw_mol, oxygen_idx)

    candidate = rw_mol.GetMol()
    Chem.SanitizeMol(candidate)
    Chem.AssignStereochemistry(candidate, force=True, cleanIt=True)
    return candidate


def _count_acetyl_groups(mol: Chem.Mol) -> int:
    patt = Chem.MolFromSmarts("[CH3][C](=O)[O]")
    return len(mol.GetSubstructMatches(patt, uniquify=True))


def _count_residual_oh(mol: Chem.Mol) -> int:
    return len(_terminal_oh_oxygen_indices(mol))


def _count_chlorine(mol: Chem.Mol) -> int:
    return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "Cl")


def build_candidates() -> list[dict]:
    seed = _load_rutinose_seed()
    seed_mol = Chem.MolFromSmiles(seed["smiles"])
    if seed_mol is None:
        raise RuntimeError("Rutinose seed SMILES could not be parsed")
    anomeric_oh_idx, anomeric_c_idx = _find_reducing_anomeric_oh(seed_mol)

    rows = []
    variants = [
        ("anomer_candidate_cw", ChiralType.CHI_TETRAHEDRAL_CW),
        ("anomer_candidate_ccw", ChiralType.CHI_TETRAHEDRAL_CCW),
    ]
    for label, tag in variants:
        candidate = _make_candidate(seed_mol, anomeric_oh_idx, anomeric_c_idx, tag)
        formula = rdMolDescriptors.CalcMolFormula(candidate)
        formula_check = "pass" if formula == EXPECTED_FORMULA else "fail"
        smiles = Chem.MolToSmiles(candidate, isomericSmiles=True)
        inchikey = rdinchi.MolToInchiKey(candidate)
        rows.append(
            {
                "candidate_id": f"reconstructed_hexa_o_acetyl_rutinosyl_chloride_{label}",
                "literature_name": TARGET_NAME,
                "reconstruction_label": label,
                "source_seed_name": seed["name"],
                "source_seed_inchikey": seed["inchikey"],
                "source_seed_smiles": seed["smiles"],
                "anomeric_chiral_tag": str(tag).replace("CHI_TETRAHEDRAL_", ""),
                "anomeric_assignment": "pending_primary_literature_validation",
                "candidate_smiles": smiles,
                "candidate_inchikey": inchikey,
                "candidate_formula": formula,
                "exact_mass": f"{Descriptors.ExactMolWt(candidate):.6f}",
                "acetyl_count": str(_count_acetyl_groups(candidate)),
                "residual_oh_count": str(_count_residual_oh(candidate)),
                "chlorine_count": str(_count_chlorine(candidate)),
                "formula_check": formula_check,
                "evidence_status": "machine_reconstructed_candidate_not_primary_confirmed",
                "stock_decision": "not_promoted_reconstruction_only",
                "template_decision": "candidate_for_mapped_template_after_primary_validation",
                "validation_notes": (
                    "Generated from true-rutinose cyclic seed by replacing reducing-end anomeric OH "
                    "with chloride and acetylating six remaining OH groups; beta assignment requires "
                    "primary-paper validation."
                ),
            }
        )
    return rows


def write_report(rows: list[dict]) -> None:
    OUTPUT_MD.parent.mkdir(parents=True, exist_ok=True)
    formula_counts = Counter(row["formula_check"] for row in rows)
    acetyl_counts = Counter(row["acetyl_count"] for row in rows)
    oh_counts = Counter(row["residual_oh_count"] for row in rows)
    chlorine_counts = Counter(row["chlorine_count"] for row in rows)
    unique_keys = {row["candidate_inchikey"] for row in rows}

    with OUTPUT_MD.open("w", encoding="utf-8") as handle:
        handle.write("# Rutinosyl Chloride Reconstruction Audit\n\n")
        handle.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        handle.write("## Summary\n\n")
        handle.write(f"- Literature target: {TARGET_NAME}\n")
        handle.write(f"- Literature source: {TARGET_DOI}\n")
        handle.write(f"- Reconstructed candidates: {len(rows)}\n")
        handle.write(f"- Unique candidate InChIKeys: {len(unique_keys)}\n")
        handle.write("- Stock promotion performed: no\n")
        handle.write("- Donor template activation performed: no\n\n")

        handle.write("## Validation Counts\n\n")
        handle.write("| Check | Value | Rows |\n")
        handle.write("|---|---|---:|\n")
        for value, count in sorted(formula_counts.items()):
            handle.write(f"| formula_check | {value} | {count} |\n")
        for value, count in sorted(acetyl_counts.items()):
            handle.write(f"| acetyl_count | {value} | {count} |\n")
        for value, count in sorted(oh_counts.items()):
            handle.write(f"| residual_oh_count | {value} | {count} |\n")
        for value, count in sorted(chlorine_counts.items()):
            handle.write(f"| chlorine_count | {value} | {count} |\n")
        handle.write("\n")

        handle.write("## Candidate Details\n\n")
        handle.write("| Candidate | Anomeric tag | Formula | InChIKey | Decision |\n")
        handle.write("|---|---|---|---|---|\n")
        for row in rows:
            handle.write(
                f"| {row['reconstruction_label']} | {row['anomeric_chiral_tag']} | "
                f"{row['candidate_formula']} | `{row['candidate_inchikey']}` | "
                f"{row['stock_decision']} |\n"
            )
        handle.write("\n")

        handle.write("## Acceptance Criteria\n\n")
        handle.write(f"- [x] Formula matches `{EXPECTED_FORMULA}` for every candidate: **{'YES' if formula_counts == {'pass': len(rows)} else 'NO'}**\n")
        handle.write(f"- [x] Six O-acetyl groups present in every candidate: **{'YES' if acetyl_counts == {'6': len(rows)} else 'NO'}**\n")
        handle.write(f"- [x] No residual OH groups remain: **{'YES' if oh_counts == {'0': len(rows)} else 'NO'}**\n")
        handle.write(f"- [x] One chloride is present: **{'YES' if chlorine_counts == {'1': len(rows)} else 'NO'}**\n")
        handle.write("- [x] Candidate anomers are not promoted to stock: **YES**\n")
        handle.write("- [x] Candidate anomers do not activate donor templates: **YES**\n\n")

        handle.write("## Interpretation\n\n")
        handle.write(
            "The reconstruction produces two anomeric machine candidates because the local true-rutinose "
            "seed does not encode a unique reducing-end anomer. These structures are appropriate for "
            "template-design and primary-literature comparison, but not for stock promotion or solved-route "
            "claims until the beta donor identity is confirmed from the primary paper or an equivalent "
            "structure source.\n"
        )


def main() -> None:
    rows = build_candidates()
    _write_csv(OUTPUT_CSV, rows, FIELDS)
    write_report(rows)
    print(f"Wrote {OUTPUT_CSV}")
    print(f"Wrote {OUTPUT_MD}")


if __name__ == "__main__":
    main()
