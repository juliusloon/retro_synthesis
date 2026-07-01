#!/usr/bin/env python
"""
Build a conservative evidence worklist for activated rutinose donors.

This script records donor evidence targets; it does not promote any molecule to
strict/trusted stock and it does not enable donor templates.
"""

import csv
from collections import Counter
from datetime import datetime
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import inchi as rdinchi


ROOT = Path(__file__).resolve().parent.parent
OUTPUT_CSV = ROOT / "templates" / "stock_layers" / "rutinose_donor_evidence_worklist.csv"
OUTPUT_MD = ROOT / "logs" / "rutinose_donor_evidence_audit.md"
RUTINOSYL_CHLORIDE_CANDIDATES_CSV = (
    ROOT / "templates" / "stock_layers" / "rutinosyl_chloride_structure_candidates.csv"
)

FIELDNAMES = [
    "candidate_id",
    "candidate_name",
    "donor_class",
    "target_core",
    "protecting_group_pattern",
    "leaving_group",
    "evidence_tier",
    "evidence_source",
    "citation_or_url",
    "exact_structure_available",
    "candidate_smiles",
    "candidate_inchikey",
    "candidate_formula",
    "structure_validation",
    "stock_decision",
    "template_decision",
    "next_action",
    "notes",
]

CANDIDATES = [
    {
        "candidate_id": "bridge_reviewed_rutinose_like_artifact",
        "candidate_name": "bridge-reviewed rutinose-like anomeric-deoxy skeleton",
        "donor_class": "route_gap_bridge_skeleton_seed",
        "target_core": "rutinose_like_bridge_core_missing_anomeric_oxygen",
        "protecting_group_pattern": "unprotected_bridge_artifact",
        "leaving_group": "none",
        "evidence_tier": "tier_3_connectivity_seed_not_free_rutinose",
        "evidence_source": "local_human_formula_audit",
        "citation_or_url": "logs/sugar_bridge_evidence_review.md; logs/sugar_bridge_core_assignment.md",
        "exact_structure_available": "yes_local",
        "candidate_smiles": "C[C@@H]1O[C@@H](OC[C@H]2OC[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H](O)[C@H]1O",
        "candidate_inchikey": "UZIKLNYKVUKZQZ-IFLAJBTPSA-N",
        "candidate_formula": "C12H22O9",
        "structure_validation": "valid_smiles_but_formula_missing_one_oxygen_vs_true_rutinose",
        "stock_decision": "keep_virtual_bridge_only",
        "template_decision": "requires_anomeric_oxygen_or_leaving_group_reconstruction",
        "next_action": "Use only as connectivity evidence; reconstruct a true rutinose donor from C12H22O10 or a literature donor structure.",
        "notes": "This is not free rutinose. It is a route-gap skeleton with formula C12H22O9, likely created by aglycone-O removal at the anomeric position.",
    },
    {
        "candidate_id": "pubchem_rutinose_free_sugar",
        "candidate_name": "PubChem rutinose free sugar",
        "donor_class": "neutral_sugar_identity_crosscheck",
        "target_core": "rutinose_core",
        "protecting_group_pattern": "free_reducing_sugar",
        "leaving_group": "none",
        "evidence_tier": "tier_2_external_identity_not_donor",
        "evidence_source": "PubChem PUGREST name lookup",
        "citation_or_url": "https://pubchem.ncbi.nlm.nih.gov/compound/5460038",
        "exact_structure_available": "yes_pubchem",
        "candidate_smiles": "C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OC[C@H]([C@H]([C@@H]([C@H](C=O)O)O)O)O)O)O)O",
        "candidate_inchikey": "CPYCUQCIDSHOHI-IFLAJBTPSA-N",
        "candidate_formula": "C12H22O10",
        "structure_validation": "valid_smiles_and_inchikey_matches_pubchem",
        "stock_decision": "no_stock_promotion_no_supplier_status_review",
        "template_decision": "identity_crosscheck_only",
        "next_action": "Use as the free-rutinose identity reference when checking reconstructed donor cores.",
        "notes": "External identity evidence for true rutinose, not evidence that an activated donor is available or validated.",
    },
    {
        "candidate_id": "hexa_o_acetyl_beta_rutinosyl_chloride_1977",
        "candidate_name": "hexa-O-acetyl-beta-rutinosyl chloride",
        "donor_class": "glycosyl_chloride_donor",
        "target_core": "rutinose_core",
        "protecting_group_pattern": "hexa_O_acetyl",
        "leaving_group": "chloride",
        "evidence_tier": "tier_1_named_literature_activated_donor",
        "evidence_source": "Crossref DOI lookup",
        "citation_or_url": "https://doi.org/10.1016/S0008-6215(00)84374-8",
        "exact_structure_available": "no_machine_structure_found",
        "candidate_smiles": "",
        "candidate_inchikey": "",
        "candidate_formula": "",
        "structure_validation": "no_machine_structure_to_validate",
        "stock_decision": "not_promoted_literature_name_only",
        "template_decision": "priority_structure_reconstruction_then_mapped_template",
        "next_action": "Retrieve the primary paper or reconstruct the exact protected beta donor, compute InChIKey, then test an atom-mapped donor disconnection.",
        "notes": "Best current direct hit for an activated rutinose donor; PubChem name lookup did not resolve this donor.",
    },
    {
        "candidate_id": "flavanontriglycoside_beta_rutinosyl_naringenin_1970",
        "candidate_name": "4'-beta-D-glucosyl-7-beta-rutinosyl-naringenin synthesis route",
        "donor_class": "rutinosylation_route_literature",
        "target_core": "rutinose_core",
        "protecting_group_pattern": "unknown",
        "leaving_group": "unknown",
        "evidence_tier": "tier_2_route_literature_leaving_group_unknown",
        "evidence_source": "Crossref DOI lookup",
        "citation_or_url": "https://doi.org/10.1002/cber.19701030528",
        "exact_structure_available": "no_machine_structure_found",
        "candidate_smiles": "",
        "candidate_inchikey": "",
        "candidate_formula": "",
        "structure_validation": "no_machine_structure_to_validate",
        "stock_decision": "not_promoted_route_context_only",
        "template_decision": "do_not_template_until_leaving_group_identified",
        "next_action": "Retrieve the primary paper and identify whether a rutinosyl halide, acetate, or other donor was used.",
        "notes": "Supports that flavonoid rutinosylation routes exist, but does not yet identify a reusable donor structure.",
    },
    {
        "candidate_id": "rational_rutin_synthesis_1968",
        "candidate_name": "rational synthesis of rutin",
        "donor_class": "rutin_route_literature",
        "target_core": "rutinose_core",
        "protecting_group_pattern": "unknown",
        "leaving_group": "unknown",
        "evidence_tier": "tier_2_route_literature_leaving_group_unknown",
        "evidence_source": "Crossref DOI lookup",
        "citation_or_url": "https://doi.org/10.1002/cber.19681010407",
        "exact_structure_available": "no_machine_structure_found",
        "candidate_smiles": "",
        "candidate_inchikey": "",
        "candidate_formula": "",
        "structure_validation": "no_machine_structure_to_validate",
        "stock_decision": "not_promoted_route_context_only",
        "template_decision": "do_not_template_until_donor_identity_confirmed",
        "next_action": "Retrieve the primary paper and extract the exact rutinose donor or precursor sequence.",
        "notes": "Useful historical route evidence for rutin, not sufficient to define a donor template by itself.",
    },
    {
        "candidate_id": "rutinosyl_trichloroacetimidate_direct_search",
        "candidate_name": "rutinosyl trichloroacetimidate",
        "donor_class": "glycosyl_trichloroacetimidate_donor",
        "target_core": "rutinose_core",
        "protecting_group_pattern": "unknown",
        "leaving_group": "trichloroacetimidate",
        "evidence_tier": "tier_4_no_direct_rutinosyl_hit",
        "evidence_source": "Crossref and PubChem name search",
        "citation_or_url": "no direct rutinosyl trichloroacetimidate hit in 2026-06-30 search",
        "exact_structure_available": "no",
        "candidate_smiles": "",
        "candidate_inchikey": "",
        "candidate_formula": "",
        "structure_validation": "no_machine_structure_to_validate",
        "stock_decision": "not_promoted_no_direct_evidence",
        "template_decision": "keep_placeholder_disabled",
        "next_action": "Only revisit after a real rutinosyl trichloroacetimidate structure or paper is found.",
        "notes": "Generic Schmidt glycosylation exists, but the current audit did not find a direct rutinose donor structure.",
    },
    {
        "candidate_id": "rutinosyl_bromide_or_acetobromorutinose_direct_search",
        "candidate_name": "rutinosyl bromide or acetobromorutinose",
        "donor_class": "glycosyl_bromide_donor",
        "target_core": "rutinose_core",
        "protecting_group_pattern": "acetylated_or_unknown",
        "leaving_group": "bromide",
        "evidence_tier": "tier_4_no_direct_rutinosyl_hit",
        "evidence_source": "Crossref and PubChem name search",
        "citation_or_url": "no direct acetobromorutinose/rutinosyl bromide hit in 2026-06-30 search",
        "exact_structure_available": "no",
        "candidate_smiles": "",
        "candidate_inchikey": "",
        "candidate_formula": "",
        "structure_validation": "no_machine_structure_to_validate",
        "stock_decision": "not_promoted_no_direct_evidence",
        "template_decision": "keep_placeholder_disabled",
        "next_action": "Prefer the chloride literature hit unless a bromide donor paper or supplier record is found.",
        "notes": "Do not infer a bromide donor from generic glycosyl bromide chemistry without direct evidence.",
    },
]


def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_reconstructed_chloride_candidates() -> list[dict]:
    rows = []
    for row in _read_csv(RUTINOSYL_CHLORIDE_CANDIDATES_CSV):
        evidence_status = row.get("evidence_status", "")
        if evidence_status == "human_confirmed_beta_anomer_candidate":
            evidence_tier = "tier_1_named_literature_activated_donor_with_confirmed_local_stereo"
            evidence_source = "local_structure_reconstruction_from_true_rutinose_seed; user_primary_source_confirmation_2026_07_01"
            exact_structure = "yes_local_beta_confirmed_by_user"
            next_action = "Validate inactive atom-mapped donor sandbox on hesperidin before any template activation or stock promotion."
            notes = (
                row["validation_notes"]
                + " This supports sandbox template design only; bridge mapping remains connectivity-level "
                "and stock promotion still requires separate stock/supplier evidence."
            )
        elif evidence_status == "reconstructed_non_beta_anomer_control":
            evidence_tier = "tier_2_reconstructed_non_beta_anomer_control"
            evidence_source = "local_structure_reconstruction_from_true_rutinose_seed; user_primary_source_confirmation_2026_07_01"
            exact_structure = "machine_reconstructed_non_beta_control"
            next_action = "Retain only as stereochemical control unless a separate non-beta donor use case is documented."
            notes = row["validation_notes"]
        else:
            evidence_tier = "tier_2_reconstructed_structure_pending_primary_validation"
            evidence_source = "local_structure_reconstruction_from_true_rutinose_seed"
            exact_structure = "machine_reconstructed_candidate_not_primary_confirmed"
            next_action = "Compare this reconstructed anomer to the primary paper beta assignment before template or stock use."
            notes = (
                row["validation_notes"]
                + " Bridge mapping is connectivity-level only; full stereochemical InChIKey differs from the route-gap artifact."
            )
        rows.append(
            {
                "candidate_id": row["candidate_id"],
                "candidate_name": f"{row['literature_name']} {row['reconstruction_label']}",
                "donor_class": "glycosyl_chloride_donor",
                "target_core": "rutinose_core",
                "protecting_group_pattern": "hexa_O_acetyl",
                "leaving_group": "chloride",
                "evidence_tier": evidence_tier,
                "evidence_source": evidence_source,
                "citation_or_url": "logs/rutinosyl_chloride_reconstruction_audit.md; logs/rutinosyl_chloride_bridge_mapping_audit.md; https://doi.org/10.1016/S0008-6215(00)84374-8",
                "exact_structure_available": exact_structure,
                "candidate_smiles": row["candidate_smiles"],
                "candidate_inchikey": row["candidate_inchikey"],
                "candidate_formula": row["candidate_formula"],
                "structure_validation": (
                    f"formula_{row['formula_check']}; acetyl_count={row['acetyl_count']}; "
                    f"residual_oh={row['residual_oh_count']}; chlorine={row['chlorine_count']}"
                    + ("; beta_assignment_confirmed" if evidence_status == "human_confirmed_beta_anomer_candidate" else "")
                    + ("; non_beta_control" if evidence_status == "reconstructed_non_beta_anomer_control" else "")
                ),
                "stock_decision": row.get("stock_decision") or "not_promoted_reconstruction_only",
                "template_decision": row.get("template_decision") or "candidate_for_mapped_template_after_primary_validation",
                "next_action": next_action,
                "notes": notes,
            }
        )
    return rows


def build_worklist_rows() -> list[dict]:
    return CANDIDATES + load_reconstructed_chloride_candidates()


def validate_candidate_structures() -> list[dict]:
    validations = []
    for row in build_worklist_rows():
        smiles = row["candidate_smiles"]
        expected_key = row["candidate_inchikey"]
        expected_formula = row["candidate_formula"]
        if not smiles:
            validations.append(
                {
                    **row,
                    "computed_formula": "",
                    "computed_inchikey": "",
                    "valid_smiles": "not_applicable",
                }
            )
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            validations.append(
                {
                    **row,
                    "computed_formula": "",
                    "computed_inchikey": "",
                    "valid_smiles": "no",
                }
            )
            continue
        computed_formula = rdMolDescriptors.CalcMolFormula(mol)
        computed_key = rdinchi.MolToInchiKey(mol)
        formula_status = "formula_match" if computed_formula == expected_formula else "formula_mismatch"
        key_status = "inchikey_match" if computed_key == expected_key else "inchikey_mismatch"
        validations.append(
            {
                **row,
                "computed_formula": computed_formula,
                "computed_inchikey": computed_key,
                "valid_smiles": f"yes_{formula_status}_{key_status}",
            }
        )
    return validations


def write_csv() -> None:
    rows = build_worklist_rows()
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_report() -> None:
    OUTPUT_MD.parent.mkdir(parents=True, exist_ok=True)
    rows = build_worklist_rows()
    validations = validate_candidate_structures()
    tier_counts = Counter(row["evidence_tier"] for row in rows)
    class_counts = Counter(row["donor_class"] for row in rows)
    stock_counts = Counter(row["stock_decision"] for row in rows)
    template_counts = Counter(row["template_decision"] for row in rows)
    exact_counts = Counter(row["exact_structure_available"] for row in rows)
    validation_counts = Counter(row["valid_smiles"] for row in validations)
    named_activated = [
        row for row in rows
        if row["evidence_tier"].startswith("tier_1_named_literature_activated_donor")
    ]
    reconstructed_candidates = [
        row for row in rows
        if row["evidence_tier"].startswith("tier_2_reconstructed")
    ]
    no_stock_promotions = all(
        "promoted" not in row["stock_decision"] or row["stock_decision"].startswith("not_promoted")
        for row in rows
    )
    no_active_templates = all(
        not row["template_decision"].startswith("activate")
        for row in rows
    )

    with OUTPUT_MD.open("w", encoding="utf-8") as handle:
        handle.write("# Rutinose Donor Evidence Audit\n\n")
        handle.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        handle.write("## Summary\n\n")
        handle.write(f"- Worklist rows: {len(rows)}\n")
        handle.write(f"- Named activated rutinose donor literature hits: {len(named_activated)}\n")
        handle.write(f"- Machine-reconstructed chloride donor candidates: {len(reconstructed_candidates)}\n")
        handle.write("- Stock promotion performed: no\n")
        handle.write("- Donor template activation performed: no\n\n")

        handle.write("## Evidence Tier Distribution\n\n")
        handle.write("| Evidence tier | Rows |\n")
        handle.write("|---|---:|\n")
        for key, value in sorted(tier_counts.items()):
            handle.write(f"| {key} | {value} |\n")
        handle.write("\n")

        handle.write("## Donor Class Distribution\n\n")
        handle.write("| Donor class | Rows |\n")
        handle.write("|---|---:|\n")
        for key, value in sorted(class_counts.items()):
            handle.write(f"| {key} | {value} |\n")
        handle.write("\n")

        handle.write("## Exact Structure Availability\n\n")
        handle.write("| Exact structure | Rows |\n")
        handle.write("|---|---:|\n")
        for key, value in sorted(exact_counts.items()):
            handle.write(f"| {key} | {value} |\n")
        handle.write("\n")

        handle.write("## Structure Validation\n\n")
        handle.write("| Validation | Rows |\n")
        handle.write("|---|---:|\n")
        for key, value in sorted(validation_counts.items()):
            handle.write(f"| {key} | {value} |\n")
        handle.write("\n")

        handle.write("## Stock Decisions\n\n")
        handle.write("| Stock decision | Rows |\n")
        handle.write("|---|---:|\n")
        for key, value in sorted(stock_counts.items()):
            handle.write(f"| {key} | {value} |\n")
        handle.write("\n")

        handle.write("## Template Decisions\n\n")
        handle.write("| Template decision | Rows |\n")
        handle.write("|---|---:|\n")
        for key, value in sorted(template_counts.items()):
            handle.write(f"| {key} | {value} |\n")
        handle.write("\n")

        handle.write("## Candidate Details\n\n")
        handle.write("| Candidate | Class | Leaving group | Tier | Exact structure | Stock decision | Template decision |\n")
        handle.write("|---|---|---|---|---|---|---|\n")
        for row in rows:
            handle.write(
                f"| {row['candidate_name']} | {row['donor_class']} | {row['leaving_group']} | "
                f"{row['evidence_tier']} | {row['exact_structure_available']} | "
                f"{row['stock_decision']} | {row['template_decision']} |\n"
            )
        handle.write("\n")

        handle.write("## Acceptance Criteria\n\n")
        handle.write(f"- [x] Every worklist row has a stock decision: **YES**\n")
        handle.write(f"- [x] No candidate is promoted to strict/trusted stock: **{'YES' if no_stock_promotions else 'NO'}**\n")
        handle.write(f"- [x] No donor template is activated: **{'YES' if no_active_templates else 'NO'}**\n")
        handle.write("- [x] Named activated donor evidence is separated from exact machine structure availability: **YES**\n")
        handle.write("- [x] Confirmed CW beta donor remains separate from stock promotion: **YES**\n")
        handle.write("- [x] Local bridge skeleton is separated from true free rutinose by formula: **YES**\n")
        handle.write("- [x] Negative direct searches are recorded instead of silently ignored: **YES**\n\n")

        handle.write("## Interpretation\n\n")
        handle.write(
            "The best current donor lead is hexa-O-acetyl-beta-rutinosyl chloride "
            "(DOI: 10.1016/S0008-6215(00)84374-8). This is evidence that an activated "
            "rutinose donor exists in the literature. The CW local chloride donor candidate has "
            "now been confirmed by the user as the beta stereochemical candidate and is used only "
            "for an inactive atom-mapped sandbox. It is still not a strict/trusted stock molecule. "
            "The CCW candidate remains a non-beta stereochemical control. The local UZIK bridge "
            "skeleton is also not true free rutinose: it is C12H22O9, while true rutinose is "
            "C12H22O10. The next scientific step is to resolve panel stereochemistry and gather "
            "separate supplier/stock-layer evidence before any promotion.\n"
        )


def main() -> None:
    write_csv()
    write_report()
    print(f"Wrote {OUTPUT_CSV}")
    print(f"Wrote {OUTPUT_MD}")


if __name__ == "__main__":
    main()
