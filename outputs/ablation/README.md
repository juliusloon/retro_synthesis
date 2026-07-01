# Ablation Outputs

This directory contains the current A1-C2 ablation outputs.

## Primary Reports

- `ablation_report.md`: Current human-readable report (A1-C2 framework).
- `ablation_report.json`: Current structured report.
- `ablation_summary.json`: Experiment run summary.
- `route_gap_worklist.csv`: Current unresolved leaf worklist.
- `route_gap_worklist.json`: JSON version of the worklist.

## Experiment Route JSONs (A1-C2)

Axis A — Discovery:
- `uspto_zinc.json` (A1)
- `uspto_rb_zinc.json` (A2)
- `uspto_custom_zinc.json` (A3)

Axis B — Evidence:
- `uspto_custom_zinc_strict.json` (B1)
- `uspto_custom_zinc_trusted.json` (B2)
- `uspto_custom_zinc_vbridge.json` (B3)

Axis C — Custom independence:
- `custom_only_zinc.json` (C1)
- `custom_only_full_stock.json` (C2)

## Legacy Outputs

Old single-target ablation outputs (baseline_*, flavonoid_*, custom_only_*) have been archived to `archive/legacy_ablation_pre_A1C2_2026-07-01/`.

## Metric Boundary

`bridge-closed` routes use `virtual_bridge` stock and are connectivity diagnostics only. Use `non-virtual effective` or strict/trusted experiments for stronger synthesis claims.
