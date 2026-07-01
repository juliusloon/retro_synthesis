# Scripts Overview

## Active Pipeline Scripts

- `custom_expansion.py`: AiZynthFinder custom expansion strategy. Applies active-template filtering, atom-ratio filtering, and post-application map-retention filtering.
- `build_reaction_family_templates.py`: Builds the active reaction-family template library from `config/reaction_families.json`.
- `run_ablation_experiments.py`: Runs the A1-C2 experiment matrix.
- `run_panel_ablation_experiments.py`: Runs the panel ablation experiments across multiple flavonoid targets.
- `compare_ablation_results.py`: Produces solved/map-valid/bridge-closed metrics and current report files.
- `compare_panel_ablation_results.py`: Produces panel-level comparison reports.
- `generate_route_gap_worklist.py`: Extracts unresolved route leaves and assigns conservative next actions.
- `build_sugar_bridge_layer.py`: Builds the route-gap-derived non-aromatic sugar bridge layer.
- `build_strict_stock.py`: Builds strict stock and separates unverified stock keys.
- `build_glycosyl_donor_surrogates.py`: Builds glycosyl donor surrogate candidates.
- `build_rutinose_donor_evidence_worklist.py`: Builds rutinose donor evidence worklist.
- `build_protected_sugar_artifact_review.py`: Builds protected sugar artifact review.
- `reconstruct_rutinosyl_chloride_candidates.py`: Reconstructs rutinosyl chloride bridge candidates.
- `validate_rutinosyl_chloride_bridge_mapping.py`: Validates rutinosyl chloride bridge mapping.
- `group_sugar_bridge_cores.py`: Groups sugar bridge core assignments.
- `fix_disaccharide_linkages.py`: Fixes disaccharide linkage analysis.
- `audit_sugar_bridge_evidence.py`: Audits sugar bridge evidence.

## Legacy Scripts (archived to `archive/legacy_scripts_2026-07-01/`)

Legacy scripts have been moved to archive. See `archive/README.md` for details.

## Current Rebuild Sequence (A1-C2 Matrix)

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_reaction_family_templates.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_ablation_experiments.py --experiments A1,A2,A3,B1,B2,B3,C1,C2
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_ablation_results.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_panel_ablation_experiments.py --experiments A1,A3,B3
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_panel_ablation_results.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/generate_route_gap_worklist.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_sugar_bridge_layer.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/audit_sugar_bridge_evidence.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/group_sugar_bridge_cores.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_protected_sugar_artifact_review.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_rutinose_donor_evidence_worklist.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_glycosyl_donor_surrogates.py
```
