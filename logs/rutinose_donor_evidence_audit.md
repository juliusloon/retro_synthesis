# Rutinose Donor Evidence Audit

Generated: 2026-07-01 16:57

## Summary

- Worklist rows: 9
- Named activated rutinose donor literature hits: 2
- Machine-reconstructed chloride donor candidates: 1
- Stock promotion performed: no
- Donor template activation performed: no

## Evidence Tier Distribution

| Evidence tier | Rows |
|---|---:|
| tier_1_named_literature_activated_donor | 1 |
| tier_1_named_literature_activated_donor_with_confirmed_local_stereo | 1 |
| tier_2_external_identity_not_donor | 1 |
| tier_2_reconstructed_non_beta_anomer_control | 1 |
| tier_2_route_literature_leaving_group_unknown | 2 |
| tier_3_connectivity_seed_not_free_rutinose | 1 |
| tier_4_no_direct_rutinosyl_hit | 2 |

## Donor Class Distribution

| Donor class | Rows |
|---|---:|
| glycosyl_bromide_donor | 1 |
| glycosyl_chloride_donor | 3 |
| glycosyl_trichloroacetimidate_donor | 1 |
| neutral_sugar_identity_crosscheck | 1 |
| route_gap_bridge_skeleton_seed | 1 |
| rutin_route_literature | 1 |
| rutinosylation_route_literature | 1 |

## Exact Structure Availability

| Exact structure | Rows |
|---|---:|
| machine_reconstructed_non_beta_control | 1 |
| no | 2 |
| no_machine_structure_found | 3 |
| yes_local | 1 |
| yes_local_beta_confirmed_by_user | 1 |
| yes_pubchem | 1 |

## Structure Validation

| Validation | Rows |
|---|---:|
| not_applicable | 5 |
| yes_formula_match_inchikey_match | 4 |

## Stock Decisions

| Stock decision | Rows |
|---|---:|
| keep_virtual_bridge_only | 1 |
| no_stock_promotion_no_supplier_status_review | 1 |
| not_promoted_literature_name_only | 1 |
| not_promoted_no_direct_evidence | 2 |
| not_promoted_no_supplier_or_stock_review | 1 |
| not_promoted_reconstruction_only | 1 |
| not_promoted_route_context_only | 2 |

## Template Decisions

| Template decision | Rows |
|---|---:|
| do_not_template_until_donor_identity_confirmed | 1 |
| do_not_template_until_leaving_group_identified | 1 |
| do_not_use_for_beta_donor_template | 1 |
| identity_crosscheck_only | 1 |
| keep_placeholder_disabled | 2 |
| priority_structure_reconstruction_then_mapped_template | 1 |
| requires_anomeric_oxygen_or_leaving_group_reconstruction | 1 |
| sandbox_mapped_template_prepared_inactive | 1 |

## Candidate Details

| Candidate | Class | Leaving group | Tier | Exact structure | Stock decision | Template decision |
|---|---|---|---|---|---|---|
| bridge-reviewed rutinose-like anomeric-deoxy skeleton | route_gap_bridge_skeleton_seed | none | tier_3_connectivity_seed_not_free_rutinose | yes_local | keep_virtual_bridge_only | requires_anomeric_oxygen_or_leaving_group_reconstruction |
| PubChem rutinose free sugar | neutral_sugar_identity_crosscheck | none | tier_2_external_identity_not_donor | yes_pubchem | no_stock_promotion_no_supplier_status_review | identity_crosscheck_only |
| hexa-O-acetyl-beta-rutinosyl chloride | glycosyl_chloride_donor | chloride | tier_1_named_literature_activated_donor | no_machine_structure_found | not_promoted_literature_name_only | priority_structure_reconstruction_then_mapped_template |
| 4'-beta-D-glucosyl-7-beta-rutinosyl-naringenin synthesis route | rutinosylation_route_literature | unknown | tier_2_route_literature_leaving_group_unknown | no_machine_structure_found | not_promoted_route_context_only | do_not_template_until_leaving_group_identified |
| rational synthesis of rutin | rutin_route_literature | unknown | tier_2_route_literature_leaving_group_unknown | no_machine_structure_found | not_promoted_route_context_only | do_not_template_until_donor_identity_confirmed |
| rutinosyl trichloroacetimidate | glycosyl_trichloroacetimidate_donor | trichloroacetimidate | tier_4_no_direct_rutinosyl_hit | no | not_promoted_no_direct_evidence | keep_placeholder_disabled |
| rutinosyl bromide or acetobromorutinose | glycosyl_bromide_donor | bromide | tier_4_no_direct_rutinosyl_hit | no | not_promoted_no_direct_evidence | keep_placeholder_disabled |
| hexa-O-acetyl-beta-rutinosyl chloride anomer_candidate_cw | glycosyl_chloride_donor | chloride | tier_1_named_literature_activated_donor_with_confirmed_local_stereo | yes_local_beta_confirmed_by_user | not_promoted_no_supplier_or_stock_review | sandbox_mapped_template_prepared_inactive |
| hexa-O-acetyl-beta-rutinosyl chloride anomer_candidate_ccw | glycosyl_chloride_donor | chloride | tier_2_reconstructed_non_beta_anomer_control | machine_reconstructed_non_beta_control | not_promoted_reconstruction_only | do_not_use_for_beta_donor_template |

## Acceptance Criteria

- [x] Every worklist row has a stock decision: **YES**
- [x] No candidate is promoted to strict/trusted stock: **YES**
- [x] No donor template is activated: **YES**
- [x] Named activated donor evidence is separated from exact machine structure availability: **YES**
- [x] Confirmed CW beta donor remains separate from stock promotion: **YES**
- [x] Local bridge skeleton is separated from true free rutinose by formula: **YES**
- [x] Negative direct searches are recorded instead of silently ignored: **YES**

## Interpretation

The best current donor lead is hexa-O-acetyl-beta-rutinosyl chloride (DOI: 10.1016/S0008-6215(00)84374-8). This is evidence that an activated rutinose donor exists in the literature. The CW local chloride donor candidate has now been confirmed by the user as the beta stereochemical candidate and is used only for an inactive atom-mapped sandbox. It is still not a strict/trusted stock molecule. The CCW candidate remains a non-beta stereochemical control. The local UZIK bridge skeleton is also not true free rutinose: it is C12H22O9, while true rutinose is C12H22O10. The next scientific step is to resolve panel stereochemistry and gather separate supplier/stock-layer evidence before any promotion.
