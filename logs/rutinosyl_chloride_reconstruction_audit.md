# Rutinosyl Chloride Reconstruction Audit

Generated: 2026-07-01 10:42

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
| anomer_candidate_cw | CW | C24H33ClO15 | `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` | not_promoted_no_supplier_or_stock_review |
| anomer_candidate_ccw | CCW | C24H33ClO15 | `BQIWKKDEVVDAAJ-KIOZKZOASA-N` | not_promoted_reconstruction_only |

## Acceptance Criteria

- [x] Formula matches `C24H33ClO15` for every candidate: **YES**
- [x] Six O-acetyl groups present in every candidate: **YES**
- [x] No residual OH groups remain: **YES**
- [x] One chloride is present: **YES**
- [x] Candidate anomers are not promoted to stock: **YES**
- [x] Candidate anomers do not activate donor templates: **YES**

## Interpretation

The reconstruction produces two anomeric machine candidates because the local true-rutinose seed does not encode a unique reducing-end anomer. The CW candidate was confirmed by the user on 2026-07-01 as the hexa-O-acetyl-beta-rutinosyl chloride stereochemical candidate, while CCW is retained only as a non-beta control. This supports inactive donor-sandbox validation, but does not promote either candidate to strict/trusted stock or activate a production template.
