# Config Overview

## Active Core Configs

- `reaction_families.json`: Source of mapped flavonoid reaction-family templates. Only entries with `active_expansion: true` are written to the active HDF5 search library.
- `sugar_closure_templates.json`: Strategy manifest for sugar closure. It includes active stock-bridge policy and disabled future template families.
- `glycosyl_donor_surrogates.json`: Glycosyl donor surrogate candidates for sugar bridge connectivity.
- `flavonoid_target_panel.csv`: Panel of 6 flavonoid targets for cross-target evaluation.
- `data.yml`: Data directory configuration.

## Active Ablation Matrix (A1-C2)

Axis A — Discovery (stock fixed to ZINC):
- `ablation_A1_uspto_zinc.yml`: USPTO only + ZINC
- `ablation_A2_uspto_rb_zinc.yml`: USPTO + RingBreaker + ZINC
- `ablation_A3_uspto_custom_zinc.yml`: USPTO + RingBreaker + Custom + ZINC

Axis B — Evidence (expansion fixed to A3):
- `ablation_B1_uspto_custom_zinc_strict.yml`: A3 + strict stock
- `ablation_B2_uspto_custom_zinc_trusted.yml`: A3 + strict + trusted stock
- `ablation_B3_uspto_custom_zinc_vbridge.yml`: A3 + strict + trusted + virtual bridge

Axis C — Custom template independence:
- `ablation_C1_custom_only_zinc.yml`: Custom only + ZINC
- `ablation_C2_custom_only_full_stock.yml`: Custom only + full stock

These are consumed by `scripts/run_ablation_experiments.py`.

## Legacy Configs (archived to `archive/legacy_configs_2026-07-01/`)

Legacy configs (baseline.yml, flavonoid.yml, hesperidin_optimized.yml) have been moved to archive. See `archive/README.md` for details.
