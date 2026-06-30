# Config Overview

## Active Core Configs

- `reaction_families.json`: Source of mapped flavonoid reaction-family templates. Only entries with `active_expansion: true` are written to the active HDF5 search library.
- `sugar_closure_templates.json`: Strategy manifest for sugar closure. It includes active stock-bridge policy and disabled future template families.

## Active Ablation Matrix

- `ablation_baseline_zinc.yml`
- `ablation_baseline_strict.yml`
- `ablation_flavonoid_zinc.yml`
- `ablation_flavonoid_strict.yml`
- `ablation_custom_only_strict.yml`
- `ablation_flavonoid_trusted.yml`
- `ablation_custom_only_trusted.yml`
- `ablation_flavonoid_virtual_bridge.yml`
- `ablation_custom_only_virtual_bridge.yml`

These are consumed by `scripts/run_ablation_experiments.py`.

## Legacy Or Compatibility Configs

- `baseline.yml`
- `flavonoid.yml`
- `hesperidin_optimized.yml`

These may be useful for older runs, but current conclusions should come from the ablation matrix and reports under `outputs/ablation/`.
