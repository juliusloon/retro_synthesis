# Scripts Overview

## Active Pipeline Scripts

- `custom_expansion.py`: AiZynthFinder custom expansion strategy. Applies active-template filtering, atom-ratio filtering, and post-application map-retention filtering.
- `build_reaction_family_templates.py`: Builds the active reaction-family template library from `config/reaction_families.json`.
- `run_ablation_experiments.py`: Runs the current nine-experiment ablation matrix.
- `compare_ablation_results.py`: Produces solved/map-valid/bridge-closed metrics and current report files.
- `generate_route_gap_worklist.py`: Extracts unresolved route leaves and assigns conservative next actions.
- `build_sugar_bridge_layer.py`: Builds the route-gap-derived non-aromatic sugar bridge layer.
- `build_strict_stock.py`: Builds strict stock and separates unverified stock keys.

## Legacy Or Support Scripts

- `run_baseline_hesperidin.py`, `run_flavonoid_hesperidin.py`, `run_optimized_hesperidin.py`: Older direct run scripts. Use only for historical reproduction.
- `compare_results.py`, `compare_all_results.py`: Older comparison helpers superseded by `compare_ablation_results.py`.
- `json2md.py`, `processing.py`: Utility/older workflow scripts.

## Current Rebuild Sequence

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_reaction_family_templates.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_ablation_experiments.py --experiments 1,2,3,4,5,6,7,8,9
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_ablation_results.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/generate_route_gap_worklist.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_sugar_bridge_layer.py --min-count 2
```
