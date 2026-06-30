# Rutinosyl Chloride Reconstruction Audit

Generated: 2026-06-30 21:29

## Summary

- Literature target: hexa-O-acetyl-beta-rutinosyl chloride
- Literature source: https://doi.org/10.1016/S0008-6215(00)84374-8
- Reconstructed candidates: 2
- Unique candidate InChIKeys: 2
- Stock promotion performed: no
- Donor template activation performed: no

## Validation Counts

| Check | Value | Rows |
|---|---|---:|
| formula_check | pass | 2 |
| acetyl_count | 6 | 2 |
| residual_oh_count | 0 | 2 |
| chlorine_count | 1 | 2 |

## Candidate Details

| Candidate | Anomeric tag | Formula | InChIKey | Decision |
|---|---|---|---|---|
| anomer_candidate_cw | CW | C24H33ClO15 | `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` | not_promoted_reconstruction_only |
| anomer_candidate_ccw | CCW | C24H33ClO15 | `BQIWKKDEVVDAAJ-KIOZKZOASA-N` | not_promoted_reconstruction_only |

## Acceptance Criteria

- [x] Formula matches `C24H33ClO15` for every candidate: **YES**
- [x] Six O-acetyl groups present in every candidate: **YES**
- [x] No residual OH groups remain: **YES**
- [x] One chloride is present: **YES**
- [x] Candidate anomers are not promoted to stock: **YES**
- [x] Candidate anomers do not activate donor templates: **YES**

## Interpretation

The reconstruction produces two anomeric machine candidates because the local true-rutinose seed does not encode a unique reducing-end anomer. These structures are appropriate for template-design and primary-literature comparison, but not for stock promotion or solved-route claims until the beta donor identity is confirmed from the primary paper or an equivalent structure source.
