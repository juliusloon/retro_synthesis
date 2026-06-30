#!/usr/bin/env python
"""
Phase 2: Fix rutinose vs neohesperidose linkage semantics.

Generates explicit structure records for:
  - Rutinose: 6-O-alpha-L-rhamnosyl-D-glucose
  - Neohesperidose: 2-O-alpha-L-rhamnosyl-D-glucose

Verifies canonical SMILES and InChIKeys differ when linkage differs.
"""

import sys
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
except ImportError:
    print("ERROR: RDKit is required.")
    sys.exit(1)


def make_rutinose_smiles():
    """
    Rutinose = 6-O-alpha-L-rhamnosyl-D-glucose
    Alpha-L-rhamnose is linked via alpha(1->6) to D-glucose.

    SMILES construction:
    - D-glucose core: OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O
    - Alpha-L-rhamnose: [C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O
    - Linkage at C6 of glucose (the CH2OH group)

    Proper rutinose SMILES with 6-O linkage:
    """
    # Rutinose: alpha-L-Rha(1->6)-D-Glc
    # The rhamnose is linked alpha to the 6-OH of glucose
    smiles = "C[C@@H]1O[C@@H](OC[C@H]2OC(O)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
    return smiles


def make_neohesperidose_smiles():
    """
    Neohesperidose = 2-O-alpha-L-rhamnosyl-D-glucose
    Alpha-L-rhamnose is linked via alpha(1->2) to D-glucose.

    SMILES construction:
    - D-glucose core with substitution at C2
    - Alpha-L-rhamnose: [C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O
    - Linkage at C2 of glucose

    Proper neohesperidose SMILES with 2-O linkage:
    """
    # Neohesperidose: alpha-L-Rha(1->2)-D-Glc
    # The rhamnose is linked alpha to the 2-OH of glucose
    smiles = "C[C@@H]1O[C@@H](O[C@H]2[C@@H](O)[C@H](O)[C@@H](O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@H]1O"
    return smiles


def verify_structure(name, smiles):
    """Verify a SMILES string is valid and compute properties."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  ERROR: Invalid SMILES for {name}: {smiles}")
        return None

    canonical = Chem.MolToSmiles(mol, canonical=True)
    inchikey = Chem.inchi.MolToInchiKey(mol) if hasattr(Chem, 'inchi') else ""

    # Count atoms
    n_atoms = mol.GetNumAtoms()
    n_rings = mol.GetRingInfo().NumRings()

    # Count hydroxyl groups
    oh_smarts = Chem.MolFromSmarts("[OX2H]")
    n_oh = len(mol.GetSubstructMatches(oh_smarts)) if oh_smarts else 0

    return {
        'name': name,
        'input_smiles': smiles,
        'canonical_smiles': canonical,
        'inchikey': inchikey,
        'n_atoms': n_atoms,
        'n_rings': n_rings,
        'n_oh': n_oh,
    }


def main():
    print("=== Phase 2: Rutinose vs Neohesperidose Linkage Fix ===\n")

    # Generate structures
    rutinose_smiles = make_rutinose_smiles()
    neohesperidose_smiles = make_neohesperidose_smiles()

    print("Generated SMILES:")
    print(f"  Rutinose:       {rutinose_smiles}")
    print(f"  Neohesperidose: {neohesperidose_smiles}")
    print()

    # Verify structures
    rutinose_info = verify_structure("rutinose", rutinose_smiles)
    neohesperidose_info = verify_structure("neohesperidose", neohesperidose_smiles)

    if not rutinose_info or not neohesperidose_info:
        print("ERROR: Could not verify structures")
        sys.exit(1)

    print("Structure verification:")
    print(f"  Rutinose:")
    print(f"    Canonical: {rutinose_info['canonical_smiles']}")
    print(f"    InChIKey:  {rutinose_info['inchikey']}")
    print(f"    Atoms: {rutinose_info['n_atoms']}, Rings: {rutinose_info['n_rings']}, OH: {rutinose_info['n_oh']}")
    print()
    print(f"  Neohesperidose:")
    print(f"    Canonical: {neohesperidose_info['canonical_smiles']}")
    print(f"    InChIKey:  {neohesperidose_info['inchikey']}")
    print(f"    Atoms: {neohesperidose_info['n_atoms']}, Rings: {neohesperidose_info['n_rings']}, OH: {neohesperidose_info['n_oh']}")
    print()

    # Check if they differ
    if rutinose_info['inchikey'] == neohesperidose_info['inchikey']:
        print("WARNING: InChIKeys are IDENTICAL - linkage not distinguished!")
        print("  This means RDKit canonicalization cannot distinguish 1->2 vs 1->6 linkage")
        print("  The generic SMILES will be used with explicit linkage annotations.")
        linkage_distinguished = False
    else:
        print("SUCCESS: InChIKeys are DIFFERENT - linkages properly distinguished!")
        linkage_distinguished = True

    print()

    # Output recommendations
    print("=== Recommendations ===\n")
    if linkage_distinguished:
        print("1. Update stock_layers_metadata.csv with distinct entries:")
        print(f"   - rutinose: InChIKey={rutinose_info['inchikey']}")
        print(f"   - neohesperidose: InChIKey={neohesperidose_info['inchikey']}")
        print()
        print("2. Update virtual_bridge_stock_inchikeys.txt with both InChIKeys")
    else:
        print("1. Keep generic SMILES but add explicit linkage annotations:")
        print(f"   - rutinose: linkage=6-O-alpha-L-rhamnosyl, confidence=medium")
        print(f"   - neohesperidose: linkage=2-O-alpha-L-rhamnosyl, confidence=medium")
        print()
        print("2. Use the same InChIKey but with different 'linkage_candidate' field")
        print("3. Mark linkage_confidence=low for both entries")

    # Save results
    output_path = Path(__file__).resolve().parent.parent / "logs" / "disaccharide_linkage_analysis.md"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("# Disaccharide Linkage Analysis\n\n")
        f.write("Generated: 2026-06-30\n\n")
        f.write("## Rutinose (6-O-alpha-L-rhamnosyl-D-glucose)\n\n")
        f.write(f"- Input SMILES: `{rutinose_smiles}`\n")
        f.write(f"- Canonical SMILES: `{rutinose_info['canonical_smiles']}`\n")
        f.write(f"- InChIKey: `{rutinose_info['inchikey']}`\n")
        f.write(f"- Atoms: {rutinose_info['n_atoms']}, Rings: {rutinose_info['n_rings']}, OH: {rutinose_info['n_oh']}\n")
        f.write(f"- Linkage: 1->6 (alpha-L-rhamnosyl to C6 of D-glucose)\n\n")

        f.write("## Neohesperidose (2-O-alpha-L-rhamnosyl-D-glucose)\n\n")
        f.write(f"- Input SMILES: `{neohesperidose_smiles}`\n")
        f.write(f"- Canonical SMILES: `{neohesperidose_info['canonical_smiles']}`\n")
        f.write(f"- InChIKey: `{neohesperidose_info['inchikey']}`\n")
        f.write(f"- Atoms: {neohesperidose_info['n_atoms']}, Rings: {neohesperidose_info['n_rings']}, OH: {neohesperidose_info['n_oh']}\n")
        f.write(f"- Linkage: 1->2 (alpha-L-rhamnosyl to C2 of D-glucose)\n\n")

        f.write("## Comparison\n\n")
        f.write(f"- Linkages distinguished by InChIKey: **{'YES' if linkage_distinguished else 'NO'}**\n")
        if linkage_distinguished:
            f.write(f"- Rutinose InChIKey: `{rutinose_info['inchikey']}`\n")
            f.write(f"- Neohesperidose InChIKey: `{neohesperidose_info['inchikey']}`\n")
        else:
            f.write("- Both structures have the same InChIKey after canonicalization\n")
            f.write("- Must use explicit linkage annotations to distinguish them\n")
        f.write("\n")

        f.write("## Bridge Layer Impact\n\n")
        if linkage_distinguished:
            f.write("- Bridge entries can be updated with distinct InChIKeys\n")
            f.write("- Rutinose-like and neohesperidose-like records will no longer share identity\n")
        else:
            f.write("- Bridge entries must use 'linkage_candidate' field to distinguish\n")
            f.write("- Both entries remain virtual_bridge with low linkage confidence\n")

    print(f"Analysis saved to: {output_path}")

    return linkage_distinguished


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
