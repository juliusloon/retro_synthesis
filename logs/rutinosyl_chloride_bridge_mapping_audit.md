# Rutinosyl Chloride Bridge Mapping Audit

Generated: 2026-06-30 21:32

## Summary

- Candidates checked: 2
- Stock promotion performed: no
- Donor template activation performed: no

## Mapping Counts

| Check | Value | Rows |
|---|---|---:|
| formula_match | yes | 2 |
| connectivity_key_match | yes | 2 |
| full_inchikey_match | no | 2 |
| mapping_decision | connectivity_match_stereochemistry_differs | 2 |

## Candidate Details

| Candidate | Collapsed formula | Collapsed InChIKey | Target InChIKey | Decision |
|---|---|---|---|---|
| reconstructed_hexa_o_acetyl_rutinosyl_chloride_anomer_candidate_cw | C12H22O9 | `UZIKLNYKVUKZQZ-WXYOOBJASA-N` | `UZIKLNYKVUKZQZ-IFLAJBTPSA-N` | connectivity_match_stereochemistry_differs |
| reconstructed_hexa_o_acetyl_rutinosyl_chloride_anomer_candidate_ccw | C12H22O9 | `UZIKLNYKVUKZQZ-WXYOOBJASA-N` | `UZIKLNYKVUKZQZ-IFLAJBTPSA-N` | connectivity_match_stereochemistry_differs |

## Acceptance Criteria

- [x] Every candidate collapses to formula `C12H22O9`: **YES**
- [x] Every candidate matches the reviewed bridge connectivity block: **YES**
- [x] Full stereochemical InChIKey match is reported separately: **YES**
- [x] No stock promotion or template activation is performed: **YES**

## Interpretation

The reconstructed chloride donor candidates collapse to the same C12H22O9 connectivity block as the reviewed route-gap bridge artifact after deacetylation and chloride-to-hydrogen normalization. Their full stereochemical InChIKeys do not match the route-gap artifact, so this is a connectivity mapping, not proof of exact bridge stereochemical identity.
