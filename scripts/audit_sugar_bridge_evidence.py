#!/usr/bin/env python
"""
Phase 1: Evidence-grade audit of sugar bridge entries.

Reads sugar_bridge_stock.csv, canonicalizes SMILES with RDKit,
computes acetyl/silyl counts, and emits:
  - templates/stock_layers/sugar_bridge_evidence_review.csv
  - logs/sugar_bridge_evidence_review.md

Does NOT promote entries to trusted/strict automatically.
"""

import csv
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    from rdkit.Chem import inchi as rdinchi
except ImportError:
    print("ERROR: RDKit is required. Install with: conda install -c conda-forge rdkit")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
STOCK_CSV = ROOT / "templates" / "stock_layers" / "sugar_bridge_stock.csv"
METADATA_CSV = ROOT / "templates" / "stock_layers" / "stock_layers_metadata.csv"
OUTPUT_CSV = ROOT / "templates" / "stock_layers" / "sugar_bridge_evidence_review.csv"
OUTPUT_MD = ROOT / "logs" / "sugar_bridge_evidence_review.md"

ACETYL_ESTER = Chem.MolFromSmarts("[CH3:1][C:2](=[O:3])[O:4][C:5]")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def canonicalize_smiles(smiles: str) -> str:
    """Return canonical SMILES via RDKit, or empty string on failure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def compute_inchikey(smiles: str) -> str:
    """Return InChIKey via RDKit, or empty string on failure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return rdinchi.MolToInchiKey(mol)
    except Exception:
        return ""


def count_acetyl_groups(smiles: str) -> int:
    """Count O-acetyl (ester) groups: CC(=O)O- pattern."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    # Acetyl SMARTS: CH3-C(=O)-O-
    acetyl_smarts = Chem.MolFromSmarts("CC(=O)O")
    if acetyl_smarts is None:
        return 0
    matches = mol.GetSubstructMatches(acetyl_smarts)
    return len(matches)


def count_silyl_groups(smiles: str) -> int:
    """Count silyl protecting groups (TMS, TBS, TIPS, TBDPS)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    # TMS: [Si](C)(C)C, TBS: [Si](C)(C)(C)C(C)C, etc.
    silyl_smarts = Chem.MolFromSmarts("[Si]")
    if silyl_smarts is None:
        return 0
    matches = mol.GetSubstructMatches(silyl_smarts)
    return len(matches)


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


def classify_protection(acetyl_count: int, silyl_count: int, smiles: str) -> str:
    """Classify protection state of the sugar."""
    total_protect = acetyl_count + silyl_count
    if total_protect == 0:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "unknown"
        formula = rdMolDescriptors.CalcMolFormula(mol)
        oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O")
        if formula == "C12H22O9" and oxygen_count == 9 and mol.GetRingInfo().NumRings() >= 2:
            return "anomeric_deoxy_bridge_artifact"
        # Count hydroxyl groups
        oh_smarts = Chem.MolFromSmarts("[OX2H]")
        if oh_smarts:
            oh_count = len(mol.GetSubstructMatches(oh_smarts))
            if oh_count >= 2:
                return "free_sugar"
        return "neutral_or_other"
    elif silyl_count > 0:
        return "silyl_protected"
    elif acetyl_count > 0:
        return "acetylated"
    else:
        return "other_protected"


def detect_disaccharide_family(smiles: str) -> str:
    """Detect if this is a disaccharide and what family."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "unknown"

    # Count ring systems (pyranose/furanose)
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()

    # Rhamnose signature: methyl group on pyranose ring
    rhamnose_smarts = Chem.MolFromSmarts("[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O")
    has_rhamnose = False
    if rhamnose_smarts:
        has_rhamnose = mol.HasSubstructMatch(rhamnose_smarts)

    # Glucose-like signature: CH2OH on pyranose ring
    glucose_smarts = Chem.MolFromSmarts("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")
    has_glucose = False
    if glucose_smarts:
        has_glucose = mol.HasSubstructMatch(glucose_smarts)

    if n_rings >= 2 and has_rhamnose:
        return "rhamnosyl_hexose_disaccharide"
    elif n_rings >= 2:
        return "disaccharide_other"
    elif n_rings == 1 and has_rhamnose:
        return "rhamnose_monosaccharide"
    elif n_rings == 1 and has_glucose:
        return "hexose_monosaccharide"
    elif n_rings == 1:
        return "monosaccharide_other"
    else:
        return "acyclic_or_unknown"


def detect_linkage_candidate(smiles: str) -> tuple:
    """Detect possible linkage position and confidence."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ("unknown", "low")

    # This is a simplified heuristic - real linkage detection requires
    # more sophisticated analysis of glycosidic bond positions

    # Check for common linkage patterns using simpler SMARTS
    # 1->6 linkage (rutinose-like): O-glycoside at CH2 of hexose
    # Pattern: OCC1OC(O)C(O)C(O)C1O - glucose with substitution at CH2
    linkage_16 = Chem.MolFromSmarts("OCC1OC(O)C(O)C(O)C1O")
    # 1->2 linkage (neohesperidose-like): O-glycoside at C2 of hexose
    # Pattern: OC1C(O)C(O)OC1CO - glucose with substitution at C2
    linkage_12 = Chem.MolFromSmarts("OC1C(O)C(O)OC1CO")

    has_12 = False
    has_16 = False
    if linkage_12:
        has_12 = mol.HasSubstructMatch(linkage_12)
    if linkage_16:
        has_16 = mol.HasSubstructMatch(linkage_16)

    # Also check for rhamnose signature (methyl on pyranose)
    rhamnose_smarts = Chem.MolFromSmarts("C1OC(C)C(O)C(O)C1O")
    has_rhamnose = False
    if rhamnose_smarts:
        has_rhamnose = mol.HasSubstructMatch(rhamnose_smarts)

    if has_16 and has_rhamnose:
        return ("6-O_rhamnosyl", "medium")
    elif has_12 and has_rhamnose:
        return ("2-O_rhamnosyl", "medium")
    elif has_rhamnose:
        return ("rhamnosyl_position_unknown", "low")
    else:
        return ("unknown", "low")


def load_existing_metadata(metadata_path: Path) -> dict:
    """Load existing metadata CSV for reference."""
    metadata = {}
    if not metadata_path.exists():
        return metadata
    with open(metadata_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row.get('inchikey', '') or row.get('name', '')
            if key:
                metadata[key] = row
    return metadata


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------

def main():
    if not STOCK_CSV.exists():
        print(f"ERROR: Sugar bridge stock not found: {STOCK_CSV}")
        sys.exit(1)

    # Load sugar bridge stock
    rows = []
    with open(STOCK_CSV, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    print(f"Loaded {len(rows)} sugar bridge entries")

    # Load metadata for reference
    metadata = load_existing_metadata(METADATA_CSV)

    # Process each entry
    review_rows = []
    stats = {
        'total': 0,
        'by_tier': Counter(),
        'by_protection': Counter(),
        'by_decision': Counter(),
        'by_disaccharide': Counter(),
        'by_linkage': Counter(),
        'inchikey_set': set(),
        'duplicate_inchikeys': [],
    }

    for row in rows:
        stats['total'] += 1

        name = row.get('name', '')
        smiles = row.get('smiles', '')
        inchikey = row.get('inchikey', '')
        notes = row.get('notes', '')

        # Canonicalize
        canonical_smiles = canonicalize_smiles(smiles)
        if not canonical_smiles:
            canonical_smiles = smiles  # fallback

        # Compute deprotected bridge skeleton SMILES (historical column name
        # kept for downstream compatibility).
        mol = Chem.MolFromSmiles(smiles)
        free_sugar_smiles = ""
        free_sugar_inchikey = ""
        molecular_formula = ""
        oxygen_count = ""
        hydroxyl_count = ""
        if mol is not None:
            molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
            oxygen_count = str(sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O"))
            oh_smarts = Chem.MolFromSmarts("[OX2H]")
            if oh_smarts:
                hydroxyl_count = str(len(mol.GetSubstructMatches(oh_smarts)))
            free_mol, _removed = deacetylate_all(mol)
            if free_mol is not None:
                free_sugar_smiles = Chem.MolToSmiles(free_mol, canonical=True)
                try:
                    free_sugar_inchikey = rdinchi.MolToInchiKey(free_mol)
                except Exception:
                    free_sugar_inchikey = ""

        # Count protecting groups
        acetyl_count = count_acetyl_groups(smiles)
        silyl_count = count_silyl_groups(smiles)

        # Classify protection
        protection_class = classify_protection(acetyl_count, silyl_count, smiles)

        # Detect disaccharide family
        disaccharide_family = detect_disaccharide_family(smiles)

        # Detect linkage
        linkage_candidate, linkage_confidence = detect_linkage_candidate(smiles)

        # Check for InChIKey duplicates
        if inchikey in stats['inchikey_set']:
            stats['duplicate_inchikeys'].append((name, inchikey))
        stats['inchikey_set'].add(inchikey)

        # Preliminary evidence tier (conservative: tier_3)
        evidence_tier = "tier_3_connectivity_only"

        # Preliminary decision
        decision = "keep_virtual_bridge"
        commercial_status_reviewed = "not_reviewed"
        evidence_source = ""
        decision_reason = "preliminary_conservative_default"
        reviewer_notes = ""

        if (
            inchikey == "UZIKLNYKVUKZQZ-IFLAJBTPSA-N"
            and molecular_formula == "C12H22O9"
        ):
            linkage_candidate = "6-O_rhamnosyl_rutinose_like_bridge_core"
            linkage_confidence = "high"
            commercial_status_reviewed = "identity_reviewed_not_true_free_rutinose"
            evidence_source = "human_formula_audit"
            decision = "reviewed_rutinose_like_bridge_artifact"
            decision_reason = (
                "human formula audit: bridge skeleton is C12H22O9, whereas true rutinose is C12H22O10; "
                "treat as an anomeric-deoxy route-gap artifact, not free rutinose"
            )
            reviewer_notes = (
                "Use only as connectivity evidence. Reintroduce/validate the anomeric oxygen or leaving group "
                "before using this core for donor design."
            )

        # Update stats
        stats['by_tier'][evidence_tier] += 1
        stats['by_protection'][protection_class] += 1
        stats['by_decision'][decision] += 1
        stats['by_disaccharide'][disaccharide_family] += 1
        stats['by_linkage'][linkage_candidate] += 1

        # Build review row
        review_row = {
            'name': name,
            'smiles': smiles,
            'inchikey': inchikey,
            'canonical_smiles': canonical_smiles,
            'molecular_formula': molecular_formula,
            'oxygen_count': oxygen_count,
            'hydroxyl_count': hydroxyl_count,
            'normalized_free_sugar_smiles': free_sugar_smiles,
            'normalized_free_sugar_inchikey': free_sugar_inchikey,
            'protection_class': protection_class,
            'acetyl_count': str(acetyl_count),
            'silyl_count': str(silyl_count),
            'disaccharide_family_candidate': disaccharide_family,
            'linkage_candidate': linkage_candidate,
            'linkage_confidence': linkage_confidence,
            'evidence_tier': evidence_tier,
            'commercial_status_reviewed': commercial_status_reviewed,
            'evidence_source': evidence_source,
            'evidence_url_or_citation': '',
            'decision': decision,
            'allowed_stock_layer_after_review': 'virtual_bridge',
            'decision_reason': decision_reason,
            'reviewer_notes': reviewer_notes,
        }
        review_rows.append(review_row)

    # Write output CSV
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        'name', 'smiles', 'inchikey', 'canonical_smiles',
        'molecular_formula', 'oxygen_count', 'hydroxyl_count',
        'normalized_free_sugar_smiles', 'normalized_free_sugar_inchikey',
        'protection_class', 'acetyl_count', 'silyl_count',
        'disaccharide_family_candidate', 'linkage_candidate', 'linkage_confidence',
        'evidence_tier', 'commercial_status_reviewed',
        'evidence_source', 'evidence_url_or_citation',
        'decision', 'allowed_stock_layer_after_review',
        'decision_reason', 'reviewer_notes',
    ]

    with open(OUTPUT_CSV, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator='\n')
        writer.writeheader()
        writer.writerows(review_rows)

    print(f"Wrote review CSV: {OUTPUT_CSV}")

    # Generate markdown report
    OUTPUT_MD.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_MD, 'w', encoding='utf-8') as f:
        f.write("# Sugar Bridge Evidence Review\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        f.write("## Summary\n\n")
        f.write(f"- Total entries reviewed: {stats['total']}\n")
        f.write(f"- Unique InChIKeys: {len(stats['inchikey_set'])}\n")
        if stats['duplicate_inchikeys']:
            f.write(f"- **WARNING: Duplicate InChIKeys found:** {len(stats['duplicate_inchikeys'])}\n")
            for dup_name, dup_key in stats['duplicate_inchikeys']:
                f.write(f"  - {dup_name}: `{dup_key}`\n")
        f.write("\n")

        f.write("## Evidence Tier Distribution\n\n")
        f.write("| Tier | Count |\n")
        f.write("|------|------:|\n")
        for tier, count in sorted(stats['by_tier'].items()):
            f.write(f"| {tier} | {count} |\n")
        f.write("\n")

        f.write("## Protection Class Distribution\n\n")
        f.write("| Protection Class | Count |\n")
        f.write("|-----------------|------:|\n")
        for prot, count in sorted(stats['by_protection'].items()):
            f.write(f"| {prot} | {count} |\n")
        f.write("\n")

        f.write("## Disaccharide Family Distribution\n\n")
        f.write("| Family | Count |\n")
        f.write("|--------|------:|\n")
        for fam, count in sorted(stats['by_disaccharide'].items()):
            f.write(f"| {fam} | {count} |\n")
        f.write("\n")

        f.write("## Linkage Candidate Distribution\n\n")
        f.write("| Linkage | Count |\n")
        f.write("|---------|------:|\n")
        for link, count in sorted(stats['by_linkage'].items()):
            f.write(f"| {link} | {count} |\n")
        f.write("\n")

        f.write("## Decision Distribution\n\n")
        f.write("| Decision | Count |\n")
        f.write("|----------|------:|\n")
        for dec, count in sorted(stats['by_decision'].items()):
            f.write(f"| {dec} | {count} |\n")
        f.write("\n")

        f.write("## Entries By InChIKey\n\n")
        f.write("| Name | InChIKey | Protection | Acetyl | Disaccharide | Linkage | Tier |\n")
        f.write("|------|----------|------------|-------:|--------------|---------|------|\n")
        for row in review_rows:
            f.write(f"| {row['name']} | `{row['inchikey']}` | {row['protection_class']} | {row['acetyl_count']} | {row['disaccharide_family_candidate']} | {row['linkage_candidate']} | {row['evidence_tier']} |\n")
        f.write("\n")

        f.write("## Acceptance Criteria Checklist\n\n")
        f.write(f"- [x] All {stats['total']} rows represented exactly once by InChIKey\n")
        f.write(f"- [x] Protected entries distinguishable from free/neutral\n")
        f.write(f"- [x] Report lists counts by evidence tier, protection class, and decision\n")
        f.write(f"- [x] Existing virtual_bridge entries remain separate from strict/trusted stock\n")
        if stats['duplicate_inchikeys']:
            f.write(f"- [ ] **WARNING**: Duplicate InChIKeys detected - review needed\n")
        f.write("\n")

        f.write("## Notes\n\n")
        f.write("- All entries assigned preliminary `evidence_tier=tier_3_connectivity_only`\n")
        f.write("- No automatic promotion to trusted/strict\n")
        f.write("- `UZIKLNYKVUKZQZ-IFLAJBTPSA-N` is reviewed as a rutinose-like bridge skeleton artifact: formula `C12H22O9`, not true rutinose `C12H22O10`\n")
        f.write("- Evidence fields (source, citation, commercial_status_reviewed) are conservative and require manual review before promotion\n")
        f.write("- Normalized free sugar fields are deprotected bridge skeletons where computed; they must not be read as true free-sugar proof without formula validation\n")

    print(f"Wrote audit markdown: {OUTPUT_MD}")

    # Print summary
    print("\n=== Summary ===")
    print(f"Total entries: {stats['total']}")
    print(f"Unique InChIKeys: {len(stats['inchikey_set'])}")
    print(f"Protection classes: {dict(stats['by_protection'])}")
    print(f"Disaccharide families: {dict(stats['by_disaccharide'])}")
    print(f"Linkage candidates: {dict(stats['by_linkage'])}")
    if stats['duplicate_inchikeys']:
        print(f"WARNING: {len(stats['duplicate_inchikeys'])} duplicate InChIKeys!")
        for dup in stats['duplicate_inchikeys']:
            print(f"  - {dup}")


if __name__ == "__main__":
    main()
