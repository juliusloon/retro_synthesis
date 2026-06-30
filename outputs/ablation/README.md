# Ablation Outputs

This directory contains the current trusted ablation outputs.

## Primary Reports

- `ablation_report.md`: Current human-readable report.
- `ablation_report.json`: Current structured report.
- `ablation_summary.json`: Experiment run summary.
- `route_gap_worklist.csv`: Current unresolved leaf worklist.
- `route_gap_worklist.json`: JSON version of the worklist.

## Experiment Route JSONs

- `baseline_zinc.json`
- `baseline_strict.json`
- `flavonoid_zinc.json`
- `flavonoid_strict.json`
- `custom_only_strict.json`
- `flavonoid_trusted.json`
- `custom_only_trusted.json`
- `flavonoid_virtual_bridge.json`
- `custom_only_virtual_bridge.json`

## Metric Boundary

`bridge-closed` routes use `virtual_bridge` stock and are connectivity diagnostics only. Use `non-virtual effective` or strict/trusted experiments for stronger synthesis claims.
