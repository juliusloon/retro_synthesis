# Sugar Bridge Evidence Review

Generated: 2026-07-01 16:56

## Summary

- Total entries reviewed: 17
- Unique InChIKeys: 17

## Evidence Tier Distribution

| Tier | Count |
|------|------:|
| tier_3_connectivity_only | 17 |

## Protection Class Distribution

| Protection Class | Count |
|-----------------|------:|
| acetylated | 16 |
| anomeric_deoxy_bridge_artifact | 1 |

## Disaccharide Family Distribution

| Family | Count |
|--------|------:|
| rhamnosyl_hexose_disaccharide | 17 |

## Linkage Candidate Distribution

| Linkage | Count |
|---------|------:|
| 6-O_rhamnosyl_rutinose_like_bridge_core | 1 |
| rhamnosyl_position_unknown | 16 |

## Decision Distribution

| Decision | Count |
|----------|------:|
| keep_virtual_bridge | 16 |
| reviewed_rutinose_like_bridge_artifact | 1 |

## Entries By InChIKey

| Name | InChIKey | Protection | Acetyl | Disaccharide | Linkage | Tier |
|------|----------|------------|-------:|--------------|---------|------|
| sugar_gap_rhamnosyl_hexose_disaccharide_UZIKLNYKVU | `UZIKLNYKVUKZQZ-IFLAJBTPSA-N` | anomeric_deoxy_bridge_artifact | 0 | rhamnosyl_hexose_disaccharide | 6-O_rhamnosyl_rutinose_like_bridge_core | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_SPCUBGINFC | `SPCUBGINFCWHNO-ZGYABWKMSA-N` | acetylated | 1 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_ASGYRTWPIY | `ASGYRTWPIYFTQG-UHBWOGNBSA-N` | acetylated | 5 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_BSESJCOUYF | `BSESJCOUYFELHO-DFRUXKTFSA-N` | acetylated | 2 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_COAOCSAWRK | `COAOCSAWRKFHQU-OYESWYITSA-N` | acetylated | 2 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_HMSBGGCYEH | `HMSBGGCYEHLLBH-NXBYORCZSA-N` | acetylated | 5 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_NDAQNSFGSP | `NDAQNSFGSPFDIX-JEUWSQIJSA-N` | acetylated | 5 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_NTQPAWBRYI | `NTQPAWBRYIOBFQ-ZGYABWKMSA-N` | acetylated | 1 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_RAFSHHVKLA | `RAFSHHVKLARFJP-YOKAVQQVSA-N` | acetylated | 5 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_SJFVSFLVJA | `SJFVSFLVJAZSFI-NXBYORCZSA-N` | acetylated | 5 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_TXZOTBUOBZ | `TXZOTBUOBZIVMS-ZGYABWKMSA-N` | acetylated | 1 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_VQFCRWDUJQ | `VQFCRWDUJQUJMR-MNENMKDJSA-N` | acetylated | 5 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_MZLKGISSOY | `MZLKGISSOYWIMT-ZFOBKIDLSA-N` | acetylated | 1 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_OFZYAYMKZL | `OFZYAYMKZLYPAR-XZTLZETGSA-N` | acetylated | 1 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_UDYLCPHGMF | `UDYLCPHGMFALPK-ZGYABWKMSA-N` | acetylated | 1 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_UNBOZMFZFL | `UNBOZMFZFLXXOQ-LIEZGGGQSA-N` | acetylated | 2 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |
| sugar_gap_rhamnosyl_hexose_disaccharide_VOAJQPZXTF | `VOAJQPZXTFIBIC-GRGQJHSESA-N` | acetylated | 2 | rhamnosyl_hexose_disaccharide | rhamnosyl_position_unknown | tier_3_connectivity_only |

## Acceptance Criteria Checklist

- [x] All 17 rows represented exactly once by InChIKey
- [x] Protected entries distinguishable from free/neutral
- [x] Report lists counts by evidence tier, protection class, and decision
- [x] Existing virtual_bridge entries remain separate from strict/trusted stock

## Notes

- All entries assigned preliminary `evidence_tier=tier_3_connectivity_only`
- No automatic promotion to trusted/strict
- `UZIKLNYKVUKZQZ-IFLAJBTPSA-N` is reviewed as a rutinose-like bridge skeleton artifact: formula `C12H22O9`, not true rutinose `C12H22O10`
- Evidence fields (source, citation, commercial_status_reviewed) are conservative and require manual review before promotion
- Normalized free sugar fields are deprotected bridge skeletons where computed; they must not be read as true free-sugar proof without formula validation
