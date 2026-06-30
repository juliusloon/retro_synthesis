#!/usr/bin/env python3
"""Build a route-gap-driven sugar bridge layer.

This script turns recurring unsolved sugar leaves into an explicit virtual
bridge stock layer. The goal is not to claim commercial availability. It gives
the search a controlled flavonoid-glycoside closure layer while keeping strict
and trusted stock semantics separate.
"""
import argparse
import csv
import re
from collections import Counter
from datetime import datetime
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, inchi, rdMolDescriptors

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ABLATION_DIR = PROJECT_ROOT / "outputs" / "ablation"
STOCK_LAYERS_DIR = PROJECT_ROOT / "templates" / "stock_layers"
LOGS_DIR = PROJECT_ROOT / "logs"

DEFAULT_GAP_WORKLIST = ABLATION_DIR / "route_gap_worklist.csv"
DEFAULT_CLUSTER_CSV = STOCK_LAYERS_DIR / "sugar_gap_clusters.csv"
DEFAULT_BRIDGE_CSV = STOCK_LAYERS_DIR / "sugar_bridge_stock.csv"
DEFAULT_METADATA = STOCK_LAYERS_DIR / "stock_layers_metadata.csv"
DEFAULT_VIRTUAL_STOCK = STOCK_LAYERS_DIR / "virtual_bridge_stock_inchikeys.txt"
DEFAULT_REPORT = LOGS_DIR / "sugar_bridge_layer_audit.md"

METADATA_FIELDS = [
    "name",
    "smiles",
    "inchikey",
    "stock_layer",
    "role",
    "evidence_level",
    "source_hint",
    "commercial_status",
    "notes",
]

CLUSTER_FIELDS = [
    "smiles",
    "canonical_smiles",
    "inchikey",
    "count",
    "first_seen_experiment",
    "upstream_family",
    "sugar_core",
    "protection_state",
    "bridge_action",
    "bridge_layer",
    "confidence",
    "mw",
    "heavy_atoms",
    "carbon_count",
    "oxygen_count",
    "sugar_ring_count",
    "aromatic_atom_count",
    "acetyl_count",
    "silyl_count",
    "rhamnose_signature",
    "notes",
]


def _mol_from_smiles(smiles: str) -> Chem.Mol | None:
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def _canonical_smiles(mol: Chem.Mol | None) -> str:
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def _inchikey(mol: Chem.Mol | None) -> str:
    if mol is None:
        return ""
    try:
        return inchi.MolToInchiKey(mol)
    except Exception:
        return ""


def _count_sugar_rings(mol: Chem.Mol) -> int:
    """Count saturated O-containing 5/6-membered C/O rings."""
    count = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6):
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(atom.GetIsAromatic() for atom in atoms):
            continue
        symbols = [atom.GetSymbol() for atom in atoms]
        if "O" not in symbols:
            continue
        if all(symbol in {"C", "O"} for symbol in symbols):
            count += 1
    return count


def _count_substructure(mol: Chem.Mol, smarts: str) -> int:
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        return 0
    return len(mol.GetSubstructMatches(patt, uniquify=True))


def _has_rhamnose_signature(mol: Chem.Mol) -> bool:
    """Heuristic for a 6-deoxy sugar methyl attached to a sugar ring."""
    patt = Chem.MolFromSmarts("[CH3]-[C;R]")
    return bool(patt is not None and mol.HasSubstructMatch(patt))


def _classify_gap(row: dict, min_count: int) -> dict:
    smiles = row.get("smiles", "")
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return {
            **row,
            "canonical_smiles": "",
            "sugar_core": "invalid",
            "protection_state": "unknown",
            "bridge_action": "reject_invalid_smiles",
            "bridge_layer": "",
            "confidence": "none",
            "mw": "",
            "heavy_atoms": "",
            "carbon_count": "",
            "oxygen_count": "",
            "sugar_ring_count": "",
            "aromatic_atom_count": "",
            "acetyl_count": "",
            "silyl_count": "",
            "rhamnose_signature": "",
            "notes": "RDKit could not parse SMILES",
        }

    canonical = _canonical_smiles(mol)
    inchikey_value = row.get("inchikey") or _inchikey(mol)
    count = int(row.get("count") or 0)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O")
    aromatic_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    sugar_ring_count = _count_sugar_rings(mol)
    acetyl_count = _count_substructure(mol, "[CH3][C](=O)[O]")
    silyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "Si")
    rhamnose_signature = _has_rhamnose_signature(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()

    if aromatic_atom_count:
        if sugar_ring_count:
            sugar_core = "flavonoid_bound_glycoside_or_aromatic_glycoside"
            bridge_action = "manual_review_do_not_bridge_as_stock"
            confidence = "medium"
            notes = "Contains aromatic scaffold; keep as chemistry target, not stock bridge"
        else:
            sugar_core = "non_sugar_aromatic"
            bridge_action = "reject_non_sugar"
            confidence = "high"
            notes = "Aromatic non-sugar gap"
    elif sugar_ring_count >= 2:
        sugar_core = (
            "rhamnosyl_hexose_disaccharide"
            if rhamnose_signature
            else "neutral_or_protected_disaccharide"
        )
        bridge_action = (
            "add_to_virtual_bridge"
            if count >= min_count
            else "manual_review_low_frequency_sugar"
        )
        confidence = "high" if count >= min_count else "medium"
        notes = "Non-aromatic disaccharide-like sugar gap"
    elif sugar_ring_count == 1:
        sugar_core = "rhamnose_like_monosaccharide" if rhamnose_signature else "hexose_like_monosaccharide"
        bridge_action = (
            "add_to_virtual_bridge"
            if count >= min_count
            else "manual_review_low_frequency_sugar"
        )
        confidence = "medium"
        notes = "Non-aromatic monosaccharide-like sugar gap"
    elif oxygen_count >= 4 and carbon_count >= 5 and not aromatic_atom_count:
        sugar_core = "acyclic_or_open_chain_sugar_like"
        bridge_action = (
            "add_to_virtual_bridge"
            if count >= min_count
            else "manual_review_low_frequency_sugar"
        )
        confidence = "medium"
        notes = "O-rich non-aromatic sugar-like gap without detected sugar ring"
    else:
        sugar_core = "non_sugar"
        bridge_action = "reject_non_sugar"
        confidence = "high"
        notes = "Does not look like a sugar bridge target"

    protection = []
    if acetyl_count:
        protection.append(f"acetylated:{acetyl_count}")
    if silyl_count:
        protection.append(f"silyl:{silyl_count}")
    if not protection:
        if sugar_core == "rhamnosyl_hexose_disaccharide" and formula == "C12H22O9":
            protection.append("anomeric_deoxy_bridge_artifact")
        else:
            protection.append("free_or_neutral")
    protection_state = ";".join(protection)

    bridge_layer = "virtual_bridge" if bridge_action == "add_to_virtual_bridge" else ""

    return {
        "smiles": smiles,
        "canonical_smiles": canonical,
        "inchikey": inchikey_value,
        "count": count,
        "first_seen_experiment": row.get("first_seen_experiment", ""),
        "upstream_family": row.get("upstream_family", ""),
        "sugar_core": sugar_core,
        "protection_state": protection_state,
        "bridge_action": bridge_action,
        "bridge_layer": bridge_layer,
        "confidence": confidence,
        "mw": f"{mw:.1f}",
        "heavy_atoms": heavy_atoms,
        "carbon_count": carbon_count,
        "oxygen_count": oxygen_count,
        "sugar_ring_count": sugar_ring_count,
        "aromatic_atom_count": aromatic_atom_count,
        "acetyl_count": acetyl_count,
        "silyl_count": silyl_count,
        "rhamnose_signature": rhamnose_signature,
        "notes": notes,
    }


def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as fileobj:
        return list(csv.DictReader(fileobj))


def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fileobj:
        writer = csv.DictWriter(fileobj, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _bridge_name(row: dict) -> str:
    core = re.sub(r"[^a-zA-Z0-9]+", "_", row["sugar_core"]).strip("_")
    suffix = (row["inchikey"] or "unknown").split("-")[0][:10]
    return f"sugar_gap_{core}_{suffix}"


def _refresh_bridge_row(row: dict) -> dict:
    """Refresh conservative labels for cumulative bridge rows."""
    refreshed = {field: row.get(field, "") for field in METADATA_FIELDS}
    mol = _mol_from_smiles(refreshed.get("smiles", ""))
    if mol is None:
        return refreshed
    formula = rdMolDescriptors.CalcMolFormula(mol)
    notes = refreshed.get("notes", "")
    if (
        formula == "C12H22O9"
        and "rhamnosyl_hexose_disaccharide" in notes
        and "free_or_neutral" in notes
    ):
        refreshed["notes"] = notes.replace(
            "free_or_neutral",
            "anomeric_deoxy_bridge_artifact",
        )
    return refreshed


def _make_bridge_rows(cluster_rows: list[dict]) -> list[dict]:
    bridge_rows = []
    seen = set()
    for row in cluster_rows:
        if row["bridge_action"] != "add_to_virtual_bridge":
            continue
        key = row["inchikey"]
        if not key or key in seen:
            continue
        seen.add(key)
        bridge_rows.append(
            {
                "name": _bridge_name(row),
                "smiles": row["canonical_smiles"] or row["smiles"],
                "inchikey": key,
                "stock_layer": "virtual_bridge",
                "role": "sugar_gap_bridge",
                "evidence_level": "route_gap_derived",
                "source_hint": "route_gap_worklist",
                "commercial_status": "not_asserted",
                "notes": (
                    f"{row['sugar_core']}; {row['protection_state']}; "
                    f"count={row['count']}; upstream={row['upstream_family']}; "
                    "virtual closure only, not strict buyable evidence"
                ),
            }
        )
    return bridge_rows


def _merge_bridge_rows(bridge_csv: Path, bridge_rows: list[dict]) -> list[dict]:
    """Merge new bridge rows with previous bridge CSV entries by InChIKey."""
    existing = _read_csv(bridge_csv)
    merged = []
    seen = set()
    # New rows come first so re-running the builder can refresh conservative
    # metadata labels without accumulating duplicate bridge entries.
    for row in bridge_rows + existing:
        key = row.get("inchikey", "")
        if not key or key in seen:
            continue
        merged.append(_refresh_bridge_row(row))
        seen.add(key)
    return merged


def _update_metadata(metadata_path: Path, bridge_rows: list[dict]) -> None:
    existing = _read_csv(metadata_path)
    new_bridge_by_key = {row["inchikey"]: row for row in bridge_rows if row.get("inchikey")}
    merged = []
    seen = set()
    for row in existing:
        key = row.get("inchikey", "")
        if key in new_bridge_by_key:
            merged.append(_refresh_bridge_row(new_bridge_by_key.pop(key)))
        else:
            merged.append(_refresh_bridge_row(row))
        if key:
            seen.add(key)
    for row in bridge_rows:
        if row["inchikey"] in seen:
            continue
        merged.append(_refresh_bridge_row(row))
        seen.add(row["inchikey"])
    _write_csv(metadata_path, merged, METADATA_FIELDS)


def _update_virtual_stock(stock_path: Path, bridge_rows: list[dict]) -> int:
    existing = []
    if stock_path.exists():
        with stock_path.open() as fileobj:
            existing = [line.strip() for line in fileobj if line.strip()]
    seen = set()
    merged = []
    for key in existing + [row["inchikey"] for row in bridge_rows]:
        if key and key not in seen:
            merged.append(key)
            seen.add(key)
    with stock_path.open("w") as fileobj:
        for key in merged:
            fileobj.write(key + "\n")
    return len(merged) - len(set(existing))


def _write_report(
    path: Path,
    gap_path: Path,
    cluster_rows: list[dict],
    bridge_rows: list[dict],
    cumulative_bridge_rows: list[dict],
    min_count: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    core_counts = Counter(row["sugar_core"] for row in cluster_rows)
    action_counts = Counter(row["bridge_action"] for row in cluster_rows)
    with path.open("w") as fileobj:
        fileobj.write("# Sugar Bridge Layer Audit\n\n")
        fileobj.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        fileobj.write(f"- Source gap worklist: `{gap_path}`\n")
        fileobj.write(f"- Min count for automatic bridge: {min_count}\n")
        fileobj.write(f"- Clustered gap leaves: {len(cluster_rows)}\n")
        fileobj.write(f"- Virtual bridge additions in this run: {len(bridge_rows)}\n")
        fileobj.write(f"- Cumulative sugar bridge entries: {len(cumulative_bridge_rows)}\n\n")

        fileobj.write("## Design Boundary\n\n")
        fileobj.write(
            "This layer only auto-bridges non-aromatic sugar-like leaves. "
            "Aromatic glycosides and flavonoid-bearing fragments are kept for "
            "manual chemistry review and are not treated as stock.\n\n"
        )

        fileobj.write("## Sugar Core Counts\n\n")
        fileobj.write("| Sugar core | Rows |\n")
        fileobj.write("|---|---:|\n")
        for core, count in core_counts.most_common():
            fileobj.write(f"| {core} | {count} |\n")
        fileobj.write("\n")

        fileobj.write("## Bridge Actions\n\n")
        fileobj.write("| Action | Rows |\n")
        fileobj.write("|---|---:|\n")
        for action, count in action_counts.most_common():
            fileobj.write(f"| {action} | {count} |\n")
        fileobj.write("\n")

        fileobj.write("## Virtual Bridge Entries Added In This Run\n\n")
        fileobj.write("| Name | InChIKey | Notes |\n")
        fileobj.write("|---|---|---|\n")
        if bridge_rows:
            rows_to_show = bridge_rows
        else:
            rows_to_show = []
        for row in rows_to_show:
            fileobj.write(f"| {row['name']} | `{row['inchikey']}` | {row['notes']} |\n")
        if not bridge_rows:
            fileobj.write("| - | - | No new bridge entries from the current worklist. |\n")

        fileobj.write("\n## Cumulative Sugar Bridge Stock\n\n")
        fileobj.write("| Name | InChIKey | Notes |\n")
        fileobj.write("|---|---|---|\n")
        for row in cumulative_bridge_rows:
            fileobj.write(f"| {row['name']} | `{row['inchikey']}` | {row['notes']} |\n")


def build_sugar_bridge_layer(
    gap_worklist: Path,
    cluster_csv: Path,
    bridge_csv: Path,
    metadata_path: Path,
    virtual_stock_path: Path,
    report_path: Path,
    min_count: int,
) -> None:
    gap_rows = _read_csv(gap_worklist)
    if not gap_rows:
        raise FileNotFoundError(f"No gap rows found in {gap_worklist}")

    cluster_rows = [_classify_gap(row, min_count=min_count) for row in gap_rows]
    cluster_rows.sort(key=lambda row: (-int(row["count"]), row["sugar_core"], row["inchikey"]))
    bridge_rows = _make_bridge_rows(cluster_rows)
    cumulative_bridge_rows = _merge_bridge_rows(bridge_csv, bridge_rows)

    _write_csv(cluster_csv, cluster_rows, CLUSTER_FIELDS)
    _write_csv(bridge_csv, cumulative_bridge_rows, METADATA_FIELDS)
    _update_metadata(metadata_path, cumulative_bridge_rows)
    added = _update_virtual_stock(virtual_stock_path, cumulative_bridge_rows)
    _write_report(
        report_path,
        gap_worklist,
        cluster_rows,
        bridge_rows,
        cumulative_bridge_rows,
        min_count,
    )

    print(f"Wrote sugar gap clusters: {cluster_csv}")
    print(f"Wrote sugar bridge stock: {bridge_csv} ({len(cumulative_bridge_rows)} cumulative rows)")
    print(f"Updated metadata: {metadata_path}")
    print(f"Updated virtual bridge stock: {virtual_stock_path} (+{added} new unique keys)")
    print(f"Wrote audit: {report_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build route-gap-driven sugar bridge stock")
    parser.add_argument("--gap-worklist", type=Path, default=DEFAULT_GAP_WORKLIST)
    parser.add_argument("--cluster-csv", type=Path, default=DEFAULT_CLUSTER_CSV)
    parser.add_argument("--bridge-csv", type=Path, default=DEFAULT_BRIDGE_CSV)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--virtual-stock", type=Path, default=DEFAULT_VIRTUAL_STOCK)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--min-count", type=int, default=2)
    args = parser.parse_args()

    build_sugar_bridge_layer(
        gap_worklist=args.gap_worklist,
        cluster_csv=args.cluster_csv,
        bridge_csv=args.bridge_csv,
        metadata_path=args.metadata,
        virtual_stock_path=args.virtual_stock,
        report_path=args.report,
        min_count=args.min_count,
    )


if __name__ == "__main__":
    main()
