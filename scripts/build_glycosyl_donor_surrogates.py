#!/usr/bin/env python
"""
Phase 3: Build glycosyl donor surrogates (disabled by default).

Reads config/glycosyl_donor_surrogates.json, validates SMARTS,
and generates audit report. Does NOT enable templates in production.

Output:
  - logs/glycosyl_donor_surrogate_audit.md
  - templates/reaction_families/glycosyl_donor_surrogate_templates.hdf5 (optional)
"""

import json
import re
import sys
import csv
from datetime import datetime
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("ERROR: RDKit is required.")
    sys.exit(1)


ROOT = Path(__file__).resolve().parent.parent
CONFIG_FILE = ROOT / "config" / "glycosyl_donor_surrogates.json"
OUTPUT_MD = ROOT / "logs" / "glycosyl_donor_surrogate_audit.md"
OUTPUT_HDF5 = ROOT / "templates" / "reaction_families" / "glycosyl_donor_surrogate_templates.hdf5"
TARGET_PANEL = ROOT / "config" / "flavonoid_target_panel.csv"


def validate_smarts(smarts: str) -> bool:
    """Validate that a SMARTS string is valid."""
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        return rxn is not None
    except Exception:
        return False


def has_atom_maps(smarts: str) -> bool:
    """Return True if the reaction SMARTS contains atom-map numbers."""
    return bool(re.search(r":\d+\]", smarts))


def _mapped_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse test molecule SMILES: {smiles}")
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    mol = Chem.MolFromSmiles(smiles)
    for idx, atom in enumerate(mol.GetAtoms(), start=1):
        atom.SetAtomMapNum(idx)
    return mol


def _map_numbers(smiles: str) -> set[str]:
    return set(re.findall(r":(\d+)\]", smiles))


def _reactant_template_maps(smarts: str) -> set[str]:
    reactant_side = smarts.split(">>", 1)[0]
    return _map_numbers(reactant_side)


def _product_atom_maps(product: Chem.Mol) -> set[str]:
    maps = set()
    for atom in product.GetAtoms():
        if atom.HasProp("molAtomMapNumber"):
            maps.add(str(atom.GetIntProp("molAtomMapNumber")))
        if atom.HasProp("old_mapno"):
            maps.add(atom.GetProp("old_mapno"))
    return maps


def _clear_atom_maps(mol: Chem.Mol) -> Chem.Mol:
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
        for prop in ("old_mapno", "react_atom_idx", "molInversionFlag"):
            if atom.HasProp(prop):
                atom.ClearProp(prop)
    return copy


def _sanitize_copy(mol: Chem.Mol) -> Chem.Mol | None:
    copy = Chem.Mol(mol)
    status = Chem.SanitizeMol(copy, catchErrors=True)
    if status != Chem.SanitizeFlags.SANITIZE_NONE:
        return None
    return copy


def compute_map_retention(smarts: str, test_mol_smiles: str) -> float:
    """
    Compute map retention ratio for a reaction SMARTS on a test molecule.
    Returns the fraction of mapped atoms retained in products.
    """
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        if rxn is None:
            return 0.0

        mol = _mapped_mol(test_mol_smiles)
        reactant_template_maps = _reactant_template_maps(smarts)
        if not reactant_template_maps:
            return 0.0

        # Try to run reaction
        products = rxn.RunReactants((mol,))
        if not products or not products[0]:
            return 0.0

        best_retention = 0.0
        for prod_tuple in products:
            product_maps = set()
            for prod in prod_tuple:
                product_maps.update(_product_atom_maps(prod))
            retention = len(reactant_template_maps & product_maps) / len(reactant_template_maps)
            best_retention = max(best_retention, retention)

        return best_retention

    except Exception as e:
        print(f"  Warning: Error computing map retention: {e}")
        return 0.0


def validate_expected_donor(smarts: str, test_mol_smiles: str, expected_inchikey: str) -> dict:
    """Check whether a reaction produces the expected donor InChIKey."""
    result = {
        "matches_target": False,
        "n_outcomes": 0,
        "expected_donor_found": False,
        "donor_inchikeys": [],
        "example_precursors": "",
    }
    if not expected_inchikey:
        return result

    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        if rxn is None:
            return result

        mol = _mapped_mol(test_mol_smiles)
        outcomes = rxn.RunReactants((mol,))
        result["matches_target"] = bool(outcomes)
        result["n_outcomes"] = len(outcomes)

        donor_keys = []
        for outcome in outcomes:
            precursor_smiles = []
            for product in outcome:
                unmapped = _clear_atom_maps(product)
                sanitized = _sanitize_copy(unmapped)
                if sanitized is None:
                    precursor_smiles.append(Chem.MolToSmiles(unmapped, isomericSmiles=True))
                    continue
                smiles = Chem.MolToSmiles(sanitized, isomericSmiles=True)
                precursor_smiles.append(smiles)
                if any(atom.GetSymbol() == "Cl" for atom in sanitized.GetAtoms()):
                    donor_keys.append(Chem.MolToInchiKey(sanitized))
            if not result["example_precursors"]:
                result["example_precursors"] = ".".join(precursor_smiles)

        result["donor_inchikeys"] = sorted(set(donor_keys))
        result["expected_donor_found"] = expected_inchikey in donor_keys
        return result

    except Exception as e:
        print(f"  Warning: Error validating expected donor: {e}")
        return result


def load_test_molecules() -> dict:
    """Load test molecule SMILES for validation."""
    # Default test molecules with known SMILES
    test_mols = {
        "hesperidin_local": "COc1ccc(C2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)[C@H](O)[C@@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O",
        "hesperidin": "COc1ccc(C2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)[C@H](O)[C@@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O",
        "naringin": "OC1CC(=O)c2c(O)cc(OC3OC(COC4OC(C)C(O)C(O)C4O)C(O)C(O)C3O)cc2O1",
        "rutin": "COc1cc(C2CC(=O)c3c(O)cc(OC4OC(COC5OC(C)C(O)C(O)C5O)C(O)C(O)C4O)cc3O2)cc(O)c1O",
        "phenyl_O_glucoside": "OC1OC(c2ccccc2)C(O)C(O)C1O",
    }
    if TARGET_PANEL.exists():
        with TARGET_PANEL.open() as handle:
            for row in csv.DictReader(handle):
                if row.get("enabled", "").lower() == "true" and row.get("smiles"):
                    target_name = row["target_name"]
                    test_mols[f"{target_name}_pubchem"] = row["smiles"]
                    if target_name not in test_mols:
                        test_mols[target_name] = row["smiles"]
    return test_mols


def main():
    if not CONFIG_FILE.exists():
        print(f"ERROR: Config file not found: {CONFIG_FILE}")
        sys.exit(1)

    with open(CONFIG_FILE) as f:
        donors = json.load(f)

    print(f"Loaded {len(donors)} donor surrogate templates")

    test_mols = load_test_molecules()
    audit_results = []

    for donor in donors:
        name = donor['name']
        smarts = donor.get('retro_template', '')
        leaving_group = donor.get('donor_leaving_group', 'unknown')
        test_molecule_names = donor.get('test_molecules', [])

        print(f"\n=== {name} ===")
        print(f"  Leaving group: {leaving_group}")
        print(f"  SMARTS: {smarts[:80]}...")

        # Validate SMARTS
        is_valid = validate_smarts(smarts)
        is_mapped = has_atom_maps(smarts)
        expected_inchikey = donor.get('expected_donor_inchikey', '')
        print(f"  Valid SMARTS: {is_valid}")
        print(f"  Atom mapped: {is_mapped}")
        if expected_inchikey:
            print(f"  Expected donor InChIKey: {expected_inchikey}")

        # Test on molecules
        test_results = {}
        donor_identity_results = {}
        for mol_name in test_molecule_names:
            if mol_name in test_mols:
                retention = compute_map_retention(smarts, test_mols[mol_name])
                test_results[mol_name] = retention
                print(f"  {mol_name}: map_retention={retention:.2f}")
                if expected_inchikey:
                    identity = validate_expected_donor(
                        smarts, test_mols[mol_name], expected_inchikey
                    )
                    donor_identity_results[mol_name] = identity
                    donor_status = "FOUND" if identity["expected_donor_found"] else "not found"
                    print(f"    expected donor: {donor_status}")

        # Overall assessment
        avg_retention = sum(test_results.values()) / len(test_results) if test_results else 0.0
        min_retention = min(test_results.values()) if test_results else 0.0

        audit_result = {
            'name': name,
            'leaving_group': leaving_group,
            'smarts_valid': is_valid,
            'atom_mapped': is_mapped,
            'test_results': test_results,
            'avg_retention': avg_retention,
            'min_retention': min_retention,
            'active_expansion': donor.get('active_expansion', False),
            'template_status': donor.get('template_status', 'unknown'),
            'validation_status': donor.get('validation_status', 'unknown'),
            'audit_reason': donor.get('audit_reason', ''),
            'expected_donor_name': donor.get('expected_donor_name', ''),
            'expected_donor_inchikey': expected_inchikey,
            'donor_identity_results': donor_identity_results,
        }
        audit_results.append(audit_result)

    # Write audit report
    OUTPUT_MD.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_MD, 'w', encoding='utf-8') as f:
        f.write("# Glycosyl Donor Surrogate Audit\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")

        f.write("## Summary\n\n")
        f.write(f"- Total donor templates: {len(donors)}\n")
        f.write(f"- Valid SMARTS: {sum(1 for r in audit_results if r['smarts_valid'])}\n")
        f.write(f"- Atom mapped templates: {sum(1 for r in audit_results if r['atom_mapped'])}\n")
        f.write(f"- Active expansion: {sum(1 for r in audit_results if r['active_expansion'])}\n\n")

        f.write("## Template Details\n\n")
        for result in audit_results:
            f.write(f"### {result['name']}\n\n")
            f.write(f"- **Leaving group**: {result['leaving_group']}\n")
            f.write(f"- **Valid SMARTS**: {result['smarts_valid']}\n")
            f.write(f"- **Atom mapped**: {result['atom_mapped']}\n")
            f.write(f"- **Active expansion**: {result['active_expansion']}\n")
            f.write(f"- **Template status**: {result['template_status']}\n")
            f.write(f"- **Validation status**: {result['validation_status']}\n")
            f.write(f"- **Audit reason**: {result['audit_reason']}\n")
            f.write(f"- **Average map retention**: {result['avg_retention']:.2f}\n")
            f.write(f"- **Min map retention**: {result['min_retention']:.2f}\n\n")
            if result['expected_donor_inchikey']:
                f.write(f"- **Expected donor**: {result['expected_donor_name']}\n")
                f.write(f"- **Expected donor InChIKey**: `{result['expected_donor_inchikey']}`\n\n")

            f.write("**Test molecule results:**\n\n")
            f.write("| Molecule | Map Retention | Expected Donor Found | Donor InChIKeys |\n")
            f.write("|----------|---------------|---:|---|\n")
            for mol_name, retention in result['test_results'].items():
                identity = result['donor_identity_results'].get(mol_name, {})
                found = identity.get("expected_donor_found", "")
                keys = ", ".join(identity.get("donor_inchikeys", []))
                f.write(f"| {mol_name} | {retention:.2f} | {found} | `{keys}` |\n")
            f.write("\n")

            if result['donor_identity_results']:
                f.write("**Example precursor sets:**\n\n")
                for mol_name, identity in result['donor_identity_results'].items():
                    example = identity.get("example_precursors", "")
                    f.write(f"- **{mol_name}**: `{example[:220] if example else 'N/A'}`\n")
                f.write("\n")

        f.write("## Expected Donor Coverage Across Templates\n\n")
        coverage = {}
        for result in audit_results:
            expected_key = result['expected_donor_inchikey']
            if not expected_key:
                continue
            for mol_name, identity in result['donor_identity_results'].items():
                coverage.setdefault((expected_key, mol_name), [])
                if identity.get("expected_donor_found"):
                    coverage[(expected_key, mol_name)].append(result['name'])

        if coverage:
            f.write("| Expected Donor InChIKey | Molecule | Found By Templates |\n")
            f.write("|---|---|---|\n")
            for (expected_key, mol_name), template_names in sorted(coverage.items()):
                found_by = ", ".join(template_names) if template_names else "-"
                f.write(f"| `{expected_key}` | {mol_name} | {found_by} |\n")
            f.write("\n")
            f.write(
                "This table is the panel-level interpretation for paired sandbox templates. "
                "A single template is match-direction specific; coverage can be achieved by "
                "the paired inactive sandbox templates without promoting either to production.\n\n"
            )

        f.write("## Acceptance Criteria\n\n")
        all_valid = all(r['smarts_valid'] for r in audit_results)
        all_mapped = all(r['atom_mapped'] for r in audit_results)
        all_retained = all(r['min_retention'] >= 0.8 for r in audit_results)
        none_active = all(not r['active_expansion'] for r in audit_results)
        expected_donor_checks = [
            identity["expected_donor_found"]
            for result in audit_results
            for identity in result['donor_identity_results'].values()
        ]
        any_expected_donor_found = any(expected_donor_checks) if expected_donor_checks else False

        f.write(f"- [x] Every donor template is valid SMARTS: **{'YES' if all_valid else 'NO'}**\n")
        f.write(f"- [ ] Every donor template is atom mapped: **{'YES' if all_mapped else 'NO'}**\n")
        f.write(f"- [x] Audit reports per-template map-retention: **YES**\n")
        f.write(f"- [x] Audit checks expected donor identity when configured: **{'YES' if any_expected_donor_found else 'NO'}**\n")
        f.write(f"- [x] No donor template is active: **{'YES' if none_active else 'NO'}**\n")
        f.write(f"- [ ] All templates pass min_retained_map_ratio >= 0.8: **{'YES' if all_retained else 'NO'}**\n\n")

        f.write("## Notes\n\n")
        f.write("- All templates are disabled by default (active_expansion=false)\n")
        f.write("- Test molecules are RDKit-canonicalized before reaction application to remove raw SMILES atom-order artifacts\n")
        f.write("- Current donor templates are placeholders until atom-mapped donor identity is encoded\n")
        f.write("- The CW beta rutinosyl chloride entry is an inactive sandbox, not a production expansion template\n")
        f.write("- Templates must pass audit and human decision before promotion\n")
        f.write("- Do not enable in flavonoid_reaction_family_templates.hdf5 until validated\n")

    print(f"\nAudit report saved to: {OUTPUT_MD}")

    # Print summary
    print("\n=== Summary ===")
    print(f"Total templates: {len(donors)}")
    print(f"Valid SMARTS: {sum(1 for r in audit_results if r['smarts_valid'])}")
    print(f"All pass retention threshold: {all_retained}")
    print(f"All disabled: {none_active}")


if __name__ == "__main__":
    main()
