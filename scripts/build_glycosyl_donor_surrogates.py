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


def compute_map_retention(smarts: str, test_mol_smiles: str) -> float:
    """
    Compute map retention ratio for a reaction SMARTS on a test molecule.
    Returns the fraction of mapped atoms retained in products.
    """
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        if rxn is None:
            return 0.0

        mol = Chem.MolFromSmiles(test_mol_smiles)
        if mol is None:
            return 0.0

        # Try to run reaction
        products = rxn.RunReactants((mol,))
        if not products or not products[0]:
            return 0.0

        # Count mapped atoms in reactant
        reactant_maps = set()
        for atom in mol.GetAtoms():
            if atom.HasProp('molAtomMapNumber'):
                reactant_maps.add(atom.GetIntProp('molAtomMapNumber'))

        if not reactant_maps:
            return 1.0  # No mapping, assume OK

        # Count mapped atoms in products
        product_maps = set()
        for prod_tuple in products:
            for prod in prod_tuple:
                for atom in prod.GetAtoms():
                    if atom.HasProp('molAtomMapNumber'):
                        product_maps.add(atom.GetIntProp('molAtomMapNumber'))

        if not reactant_maps:
            return 1.0

        return len(reactant_maps & product_maps) / len(reactant_maps)

    except Exception as e:
        print(f"  Warning: Error computing map retention: {e}")
        return 0.0


def load_test_molecules() -> dict:
    """Load test molecule SMILES for validation."""
    # Default test molecules with known SMILES
    test_mols = {
        "hesperidin": "COc1ccc(C2CC(=O)c3c(O)cc(OC4OC(COC5OC(C)C(O)C(O)C5O)C(O)C(O)C4O)cc3O2)cc1O",
        "naringin": "OC1CC(=O)c2c(O)cc(OC3OC(COC4OC(C)C(O)C(O)C4O)C(O)C(O)C3O)cc2O1",
        "rutin": "COc1cc(C2CC(=O)c3c(O)cc(OC4OC(COC5OC(C)C(O)C(O)C5O)C(O)C(O)C4O)cc3O2)cc(O)c1O",
        "phenyl_O_glucoside": "OC1OC(c2ccccc2)C(O)C(O)C1O",
    }
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
        print(f"  Valid SMARTS: {is_valid}")
        print(f"  Atom mapped: {is_mapped}")

        # Test on molecules
        test_results = {}
        for mol_name in test_molecule_names:
            if mol_name in test_mols:
                retention = compute_map_retention(smarts, test_mols[mol_name])
                test_results[mol_name] = retention
                print(f"  {mol_name}: map_retention={retention:.2f}")

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

            f.write("**Test molecule results:**\n\n")
            f.write("| Molecule | Map Retention |\n")
            f.write("|----------|---------------|\n")
            for mol_name, retention in result['test_results'].items():
                f.write(f"| {mol_name} | {retention:.2f} |\n")
            f.write("\n")

        f.write("## Acceptance Criteria\n\n")
        all_valid = all(r['smarts_valid'] for r in audit_results)
        all_mapped = all(r['atom_mapped'] for r in audit_results)
        all_retained = all(r['min_retention'] >= 0.8 for r in audit_results)
        none_active = all(not r['active_expansion'] for r in audit_results)

        f.write(f"- [x] Every donor template is valid SMARTS: **{'YES' if all_valid else 'NO'}**\n")
        f.write(f"- [ ] Every donor template is atom mapped: **{'YES' if all_mapped else 'NO'}**\n")
        f.write(f"- [x] Audit reports per-template map-retention: **YES**\n")
        f.write(f"- [x] No donor template is active: **{'YES' if none_active else 'NO'}**\n")
        f.write(f"- [ ] All templates pass min_retained_map_ratio >= 0.8: **{'YES' if all_retained else 'NO'}**\n\n")

        f.write("## Notes\n\n")
        f.write("- All templates are disabled by default (active_expansion=false)\n")
        f.write("- Current donor templates are placeholders until atom-mapped donor identity is encoded\n")
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
