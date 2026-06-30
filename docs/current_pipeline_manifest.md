# Current Retrosynthesis Pipeline Manifest

Last updated: 2026-06-30

This manifest defines the active hesperidin/flavonoid-glycoside retrosynthesis workflow. It is intended to prevent legacy template experiments, virtual bridge hypotheses, and strict stock evidence from being mixed together.

## Current Objective

Use AiZynthFinder plus audited flavonoid reaction-family templates to evaluate hesperidin-like flavonoid glycoside retrosynthesis.

The current scientific distinction is:

- **strict/trusted solved**: meaningful synthesis evidence candidate.
- **bridge-closed solved**: connectivity validation using `virtual_bridge`; useful for diagnosis, not proof of buyability.
- **ZINC solved**: useful baseline, often closes by treating complex database molecules as stock.

## Active Pipeline

Run from the repository root:

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_reaction_family_templates.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_ablation_experiments.py --experiments 1,2,3,4,5,6,7,8,9
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_ablation_results.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/generate_route_gap_worklist.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_sugar_bridge_layer.py --min-count 2
```

For quick sugar-bridge validation:

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_ablation_experiments.py --experiments 8,9
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_ablation_results.py
```

## Active Scripts

| Path | Status | Role |
|---|---|---|
| `scripts/custom_expansion.py` | active | Uniform custom AiZynthFinder expansion strategy with active-template filtering and map-retention post-filtering. |
| `scripts/build_reaction_family_templates.py` | active | Builds audited reaction-family HDF5/CSV from `config/reaction_families.json`. Only active templates are written to search libraries. |
| `scripts/run_ablation_experiments.py` | active | Runs the nine-experiment stock/template ablation matrix. |
| `scripts/compare_ablation_results.py` | active | Computes AiZynth solved, map-valid solved, bridge-closed solved, stock-layer diagnostics, and family usage. |
| `scripts/generate_route_gap_worklist.py` | active | Generates unresolved leaf worklists after ablation runs. |
| `scripts/build_sugar_bridge_layer.py` | active | Builds gap-derived non-aromatic sugar virtual bridge entries and audit tables. |
| `scripts/build_strict_stock.py` | active/support | Builds strict buyable stock and separates unverified keys. |

## Active Configs

| Path | Status | Role |
|---|---|---|
| `config/reaction_families.json` | active | Source of audited flavonoid reaction-family templates. |
| `config/sugar_closure_templates.json` | active strategy manifest | Documents sugar closure strategy and future disabled template families. |
| `config/ablation_*.yml` | active experiments | Nine ablation configurations for ZINC/strict/trusted/virtual bridge tests. |
| `config/flavonoid.yml` | legacy/compatibility | Older flavonoid configuration; do not use for current conclusions unless explicitly refreshed. |
| `config/hesperidin_optimized.yml` | legacy/compatibility | Older optimized hesperidin configuration; do not use for current conclusions unless explicitly refreshed. |
| `config/baseline.yml` | legacy baseline | Older baseline config. |

## Active Template And Stock Artifacts

| Path | Status | Role |
|---|---|---|
| `templates/reaction_families/flavonoid_reaction_family_templates.hdf5` | active generated | MCTS search library containing active reaction-family templates only. |
| `templates/reaction_families/flavonoid_reaction_family_templates_audit.csv` | active audit | Full active/inactive template audit table. |
| `templates/stock_layers/strict_buyable_stock_inchikeys.txt` | active stock | Strict stock layer. |
| `templates/stock_layers/trusted_intermediate_stock_inchikeys.txt` | active stock | Trusted intermediate layer. |
| `templates/stock_layers/virtual_bridge_stock_inchikeys.txt` | active diagnostic stock | Virtual bridge stock for connectivity diagnosis only. |
| `templates/stock_layers/sugar_gap_clusters.csv` | active audit | Classification of route-gap sugar leaves. |
| `templates/stock_layers/sugar_bridge_stock.csv` | active audit/generated stock | Route-gap-derived sugar bridge entries. |
| `templates/stock_layers/stock_layers_metadata.csv` | active metadata | Stock-layer evidence and role metadata. |

## Current Trusted Outputs

| Path | Status | Role |
|---|---|---|
| `outputs/ablation/ablation_report.md` | active report | Human-readable latest ablation summary. |
| `outputs/ablation/ablation_report.json` | active report | Structured latest ablation summary. |
| `outputs/ablation/route_gap_worklist.csv` | active worklist | Latest unresolved leaf worklist after current bridge layer. |
| `logs/reaction_family_template_audit.md` | active audit | Reaction-family template audit. |
| `logs/sugar_bridge_layer_audit.md` | active audit | Sugar bridge layer audit. |

## Legacy Areas

| Path | Status | Notes |
|---|---|---|
| `src/` | legacy template experiments | Early extraction/conversion/diagnostic scripts. Keep for provenance, but do not use for current conclusions without review. |
| `templates/flavonoid_templates.*` | legacy custom library | Older template set; current MCTS ablations use reaction-family templates. |
| `templates/flavonoid_structural_diversity/` | legacy template family | Structural diversity templates retained for provenance. |
| `templates/flavonoid_biosynthesis/` | legacy template family | Biosynthesis templates retained for provenance. |
| `templates/literature_curated/` | evidence/provenance | Literature-derived templates and evidence tables. `reaction_family_evidence_minimal.csv` is still relevant evidence. |
| `archive/` | archived outputs/logs | Older route outputs and dated reports moved out of active paths. |

## Current Result Snapshot

After sugar bridge layer construction:

| Experiment | Effective solved | Bridge-closed | Non-virtual |
|---|---:|---:|---:|
| `custom_only_virtual_bridge` | 2 | 2 | 0 |
| `flavonoid_virtual_bridge` | 10 | 10 | 0 |
| `flavonoid_strict` | 0 | 0 | 0 |
| `flavonoid_trusted` | 0 | 0 | 0 |

Interpretation: the audited glycoside cleavage template reaches the correct sugar/aglycone disconnection. The solved virtual-bridge routes validate connectivity but do not prove strict buyability.

## Next Scientific Tasks

1. Evidence-grade the 17 sugar bridge entries in `templates/stock_layers/sugar_bridge_stock.csv`.
2. Separate rutinose-like and neohesperidose-like linkages more rigorously.
3. Build and audit activated glycosyl donor surrogate templates only after map-retention validation.
4. Keep aromatic flavonoid glycoside leaves as manual review targets, not virtual stock.
5. Add a flavonoid glycoside test panel beyond hesperidin.
