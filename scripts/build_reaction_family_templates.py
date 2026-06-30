#!/usr/bin/env python3
"""Build atom-mapped flavonoid reaction-family templates.

The existing structural and biosynthesis templates are useful as domain priors,
but many are scaffold-collapse transformations. This script creates a small
reaction-family template library whose SMARTS carry atom maps and can be gated
by post-application map retention in scripts/custom_expansion.py.
"""
import argparse
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FAMILIES = PROJECT_ROOT / "config" / "reaction_families.json"
DEFAULT_OUTPUT = PROJECT_ROOT / "templates" / "reaction_families" / "flavonoid_reaction_family_templates.hdf5"
DEFAULT_REPORT = PROJECT_ROOT / "logs" / "reaction_family_template_audit.md"
DEFAULT_AUDIT_TABLE = PROJECT_ROOT / "templates" / "reaction_families" / "flavonoid_reaction_family_templates_audit.csv"

# Test molecules for validation
TEST_MOLECULES = {
    "hesperidin": "COc1ccc(C2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)[C@H](O)[C@@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O",
    "hesperetin": "COc1ccc(C2CC(=O)c3c(O)cc(O)cc3O2)cc1O",
    "naringenin": "OC1CC(=O)c2c(O)cc(O)cc2O1",
    "hydroxychalcone_2": "O=C(c1ccccc1)/C=C/c1ccccc1O",
    "simple_chalcone": "O=C(c1ccccc1)/C=C/c1ccccc1",
}


def _mapped_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse target SMILES: {smiles}")
    for idx, atom in enumerate(mol.GetAtoms(), start=1):
        atom.SetAtomMapNum(idx)
    return mol


def _map_numbers(smiles: str) -> set[str]:
    return set(re.findall(r":(\d+)\]", smiles))


def _validate_template(retro_template: str, target_smiles: str) -> dict:
    """Validate a single template against a target molecule."""
    rxn = AllChem.ReactionFromSmarts(retro_template)
    if rxn is None:
        return {
            "valid_smarts": False,
            "matches_target": False,
            "n_outcomes": 0,
            "max_map_retention": 0.0,
            "example_precursors": "",
        }

    target = _mapped_mol(target_smiles)
    target_maps = _map_numbers(Chem.MolToSmiles(target, isomericSmiles=True))
    outcomes = rxn.RunReactants((target,))
    best_retention = 0.0
    best_precursors = ""

    for outcome in outcomes:
        precursor_smiles = ".".join(
            Chem.MolToSmiles(mol, isomericSmiles=True) for mol in outcome
        )
        precursor_maps = _map_numbers(precursor_smiles)
        retention = (
            len(target_maps & precursor_maps) / len(target_maps)
            if target_maps
            else 1.0
        )
        if retention > best_retention:
            best_retention = retention
            best_precursors = precursor_smiles

    return {
        "valid_smarts": True,
        "matches_target": bool(outcomes),
        "n_outcomes": len(outcomes),
        "max_map_retention": best_retention,
        "example_precursors": best_precursors,
    }


def _validate_template_all_targets(retro_template: str, test_molecules: dict) -> dict:
    """Validate a template against all test molecules."""
    results = {}
    for name, smiles in test_molecules.items():
        results[name] = _validate_template(retro_template, smiles)
    return results


def build_templates(families_path: Path, output_path: Path, report_path: Path) -> None:
    """Build templates and audit against all test molecules."""
    with families_path.open() as fileobj:
        families = json.load(fileobj)

    rows = []
    audit_rows = []
    for family in families:
        template = family["retro_template"]
        template_hash = hashlib.sha256(template.encode("utf-8")).hexdigest()
        active_expansion = bool(family.get("active_expansion", True))

        # Validate against all test molecules
        all_audits = _validate_template_all_targets(template, TEST_MOLECULES)

        # Get best retention across all test molecules
        best_retention = max(audit["max_map_retention"] for audit in all_audits.values())
        best_outcomes = max(audit["n_outcomes"] for audit in all_audits.values())

        # Check if passes min_retained_map_ratio
        min_ratio = family.get("min_retained_map_ratio", 0.8)
        passes_min_ratio = best_retention >= min_ratio

        rows.append(
            {
                "retro_template": template,
                "template_hash": template_hash,
                "classification": family["classification"],
                "template_class": family.get("template_class", "unknown"),
                "family_name": family["name"],
                "evidence_family": family.get("evidence_family", ""),
                "notes": family.get("notes", ""),
                "active_expansion": active_expansion,
                "audit_reason": family.get("audit_reason", ""),
                "max_map_retention": best_retention,
                "max_outcomes": best_outcomes,
                "passes_min_retained_map_ratio": passes_min_ratio,
                "min_retained_map_ratio_threshold": min_ratio,
            }
        )
        audit_rows.append((family, template_hash, all_audits))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    full_df = pd.DataFrame(rows)
    active_df = full_df[full_df["active_expansion"]].reset_index(drop=True)
    active_df.to_hdf(output_path, key="table", mode="w", format="table")
    active_df.to_csv(output_path.with_suffix(".csv"), index=True, sep="\t")
    full_df.to_csv(DEFAULT_AUDIT_TABLE, index=False)

    # Generate detailed audit report
    with report_path.open("w") as fileobj:
        fileobj.write("# Reaction Family Template Audit\n\n")
        fileobj.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        fileobj.write(f"- Source: `{families_path}`\n")
        fileobj.write(f"- Output HDF5: `{output_path}`\n")
        fileobj.write(f"- Output CSV: `{output_path.with_suffix('.csv')}`\n\n")
        fileobj.write(
            f"- Audit CSV: `{DEFAULT_AUDIT_TABLE}`\n"
            f"- Families in source: {len(full_df)}\n"
            f"- Active templates written to search library: {len(active_df)}\n"
            f"- Inactive templates retained only for audit: {len(full_df) - len(active_df)}\n\n"
        )

        fileobj.write("## Test Molecules\n\n")
        for name, smiles in TEST_MOLECULES.items():
            fileobj.write(f"- **{name}**: `{smiles}`\n")
        fileobj.write("\n")

        fileobj.write("## Template Summary\n\n")
        fileobj.write("| Family | Active | Classification | Template Class | Valid SMARTS | Max Map Retention | Passes Min Ratio | Evidence Family |\n")
        fileobj.write("|---|---:|---|---|---:|---:|---:|---|\n")
        for family, template_hash, all_audits in audit_rows:
            best_retention = max(audit["max_map_retention"] for audit in all_audits.values())
            min_ratio = family.get("min_retained_map_ratio", 0.8)
            passes = "✓" if best_retention >= min_ratio else "✗"
            active = "✓" if family.get("active_expansion", True) else "✗"
            fileobj.write(
                f"| {family['name']} | {active} | {family['classification']} | "
                f"{family.get('template_class', 'unknown')} | "
                f"{all_audits['hesperidin']['valid_smarts']} | "
                f"{best_retention:.3f} | {passes} | "
                f"{family.get('evidence_family', '')} |\n"
            )
        fileobj.write("\n")

        fileobj.write("## Detailed Validation Results\n\n")
        for family, template_hash, all_audits in audit_rows:
            fileobj.write(f"### {family['name']}\n\n")
            fileobj.write(f"- **Classification**: {family['classification']}\n")
            fileobj.write(f"- **Template Class**: {family.get('template_class', 'unknown')}\n")
            fileobj.write(f"- **Evidence Family**: {family.get('evidence_family', '')}\n")
            fileobj.write(f"- **Active Expansion**: {family.get('active_expansion', True)}\n")
            fileobj.write(f"- **Audit Reason**: {family.get('audit_reason', '')}\n")
            fileobj.write(f"- **Min Retained Map Ratio**: {family.get('min_retained_map_ratio', 0.8)}\n")
            fileobj.write(f"- **Notes**: {family.get('notes', '')}\n\n")

            fileobj.write("| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |\n")
            fileobj.write("|---|---|---|---:|---:|---|\n")
            for mol_name, audit in all_audits.items():
                fileobj.write(
                    f"| {mol_name} | {audit['valid_smarts']} | "
                    f"{audit['matches_target']} | {audit['n_outcomes']} | "
                    f"{audit['max_map_retention']:.3f} | "
                    f"`{audit['example_precursors'][:50] if audit['example_precursors'] else 'N/A'}` |\n"
                )
            fileobj.write("\n")

        fileobj.write("## Notes\n\n")
        fileobj.write(
            "Only templates with `active_expansion: true` are written to the HDF5 "
            "and CSV search libraries. Inactive templates are retained in the "
            "audit report and audit CSV for traceability.\n"
        )
        fileobj.write("\n")
        fileobj.write(
            "The expansion strategy also gates applied actions with "
            "`min_retained_map_ratio`, so target-validation retention and "
            "post-application route retention are both visible.\n"
        )
        fileobj.write("\n")
        fileobj.write(
            "All templates must have valid SMARTS and atom maps. "
            "Templates with 0 map retention are flagged for review.\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Build mapped flavonoid reaction-family templates")
    parser.add_argument("--families", type=Path, default=DEFAULT_FAMILIES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    build_templates(args.families, args.output, args.report)
    print(f"Wrote templates: {args.output}")
    print(f"Wrote CSV: {args.output.with_suffix('.csv')}")
    print(f"Wrote audit CSV: {DEFAULT_AUDIT_TABLE}")
    print(f"Wrote audit: {args.report}")


if __name__ == "__main__":
    main()
