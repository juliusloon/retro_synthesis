# Next AI Execution Plan: From Virtual Sugar Bridge To Auditable Donor Layer

Last updated: 2026-06-30

This document is an executable handoff for the best AI agent. It is based on
`docs/current_pipeline_manifest.md` and the current active pipeline. The main
problem is not route connectivity anymore: `flavonoid_virtual_bridge` is
`10/10` effective solved, but every solved route is `bridge-closed`. The
scientific task is to turn the sugar bridge diagnosis into auditable chemistry
without weakening strict/trusted solved definitions.

## Non-Negotiable Boundaries

1. Do not count `virtual_bridge` leaves as strict or trusted solved evidence.
2. Do not add aromatic flavonoid glycoside leaves to stock automatically.
3. Do not use legacy `src/` scripts or legacy template folders for current
   conclusions unless explicitly refreshing and auditing them.
4. Preserve the split among:
   - `strict/trusted solved`
   - `bridge-closed solved`
   - `ZINC solved`
5. Do not enable chalcone templates until they pass map-retention validation on
   representative test molecules.
6. Do not treat protected sugar artifacts from USPTO as real donor evidence
   without normalization and provenance review.

## Current Baseline To Verify Before Editing

From the repository root, run:

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_reaction_family_templates.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_ablation_experiments.py --experiments 1,2,3,4,5,6,7,8,9
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_ablation_results.py
```

Expected current signal:

- `flavonoid_virtual_bridge`: effective solved `10`, bridge-closed `10`,
  non-virtual `0`.
- `flavonoid_strict`: effective solved `0`.
- `flavonoid_trusted`: effective solved `0`.
- `templates/stock_layers/sugar_bridge_stock.csv` has 17 bridge entries.
- `outputs/ablation/route_gap_worklist.csv` top remaining gaps are aromatic
  glycoside/manual-review leaves after sugar-core closure.

If these values differ, record the new baseline in a short note before making
scientific changes.

## Phase 1: Evidence-Grade The 17 Sugar Bridge Entries

Goal: create an auditable review table. Do not update trusted/strict stock in
this phase unless there is explicit evidence and the decision is separable.

Input:

- `templates/stock_layers/sugar_bridge_stock.csv`
- `templates/stock_layers/stock_layers_metadata.csv`
- `logs/sugar_bridge_layer_audit.md`
- `outputs/ablation/route_gap_worklist.csv`

Create:

- `templates/stock_layers/sugar_bridge_evidence_review.csv`
- `logs/sugar_bridge_evidence_review.md`

Required CSV columns:

```text
name
smiles
inchikey
canonical_smiles
normalized_free_sugar_smiles
normalized_free_sugar_inchikey
protection_class
acetyl_count
silyl_count
disaccharide_family_candidate
linkage_candidate
linkage_confidence
evidence_tier
commercial_status_reviewed
evidence_source
evidence_url_or_citation
decision
allowed_stock_layer_after_review
decision_reason
reviewer_notes
```

Use these evidence tiers:

- `tier_0_strict_commercial`: exact compound or accepted starting material is
  commercially available with a verifiable supplier/catalog or equivalent.
- `tier_1_literature_verified`: exact compound or donor/intermediate is
  reported in literature with adequate structural identity.
- `tier_2_reasonable_intermediate`: chemically plausible protected/free sugar
  intermediate, but not enough exact evidence for trusted/strict solved.
- `tier_3_connectivity_only`: useful only as virtual bridge connectivity.
- `tier_4_artifact_reject`: likely template/protection artifact; exclude from
  bridge or donor conclusions.

Use these decisions:

- `keep_virtual_bridge`
- `promote_to_trusted_intermediate_candidate`
- `promote_to_strict_candidate`
- `normalize_to_free_sugar_core`
- `reject_artifact`
- `needs_human_review`

Minimum automation:

1. Write or update a small script, preferably
   `scripts/audit_sugar_bridge_evidence.py`, that reads
   `sugar_bridge_stock.csv`, canonicalizes SMILES with RDKit, computes
   acetyl/silyl counts, and emits the review CSV with blank evidence fields.
2. The script may assign preliminary `protection_class` and conservative
   default `evidence_tier=tier_3_connectivity_only`.
3. It must not promote entries automatically to trusted/strict.

Acceptance criteria:

- The 17 rows are all represented exactly once by InChIKey.
- Protected entries are distinguishable from free/neutral entries.
- The report lists counts by evidence tier, protection class, and decision.
- Existing `virtual_bridge` entries remain separate from strict/trusted stock.

## Phase 2: Fix Rutinose vs Neohesperidose Linkage Semantics

Goal: remove the current metadata ambiguity where rutinose and neohesperidose
are represented with the same InChIKey/SMILES.

Files to inspect and update:

- `templates/stock_layers/stock_layers_metadata.csv`
- `templates/stock_layers/virtual_bridge_stock_inchikeys.txt`
- `config/sugar_closure_templates.json`
- Any generated review/audit files from Phase 1

Work:

1. Generate explicit structure records for:
   - rutinose: 6-O-alpha-L-rhamnosyl-D-glucose-like linkage
   - neohesperidose: 2-O-alpha-L-rhamnosyl-D-glucose-like linkage
2. Verify that canonical SMILES and InChIKeys differ when linkage differs.
3. If stereochemistry is uncertain, mark `linkage_confidence=low` and keep the
   row virtual only.
4. Update notes so future agents cannot infer exact linkage from the generic
   `rhamnosyl_hexose_disaccharide` class alone.

Acceptance criteria:

- Rutinose-like and neohesperidose-like records no longer share the same
  exact identity unless a deliberate generic placeholder is explicitly named as
  such.
- The bridge layer can still run without crashing.
- The report states which bridge entries can be assigned to 2-O, 6-O, generic,
  or unknown linkage classes.

## Phase 3: Build A Donor-Surrogate Sandbox, Disabled By Default

Goal: model real glycosylation chemistry as donor formation/use, but keep it
out of production MCTS until map-retention and identity audits pass.

Do not enable the existing `activated_glycosyl_donor_disconnection` template as
is. It currently has the same SMARTS as `o_glycoside_cleavage_pyranose` and
does not encode a leaving group.

Create or update:

- `config/glycosyl_donor_surrogates.json`
- `scripts/build_glycosyl_donor_surrogates.py`
- `logs/glycosyl_donor_surrogate_audit.md`
- optionally
  `templates/reaction_families/glycosyl_donor_surrogate_templates.hdf5`

Candidate donor families to encode as disabled templates first:

- glycosyl bromide donor
- glycosyl chloride donor
- glycosyl trichloroacetimidate donor
- thioglycoside donor
- glycosyl acetate or carbonate surrogate only if clearly marked as surrogate

Required template properties:

```json
{
  "name": "...",
  "classification": "glycosyl_donor_surrogate",
  "retro_template": "...",
  "donor_leaving_group": "...",
  "target_linkage_scope": "aryl_O_glycoside|disaccharide|generic",
  "allowed_stock_layers": ["trusted_intermediate", "virtual_bridge"],
  "min_retained_map_ratio": 0.8,
  "active_expansion": false,
  "audit_reason": "disabled until mapped donor identity and test molecules pass"
}
```

Required test molecules:

- hesperidin
- naringin or neohesperidin-like glycoside
- rutin-like glycoside
- one simple phenyl O-glycoside
- one donor-only test molecule for each leaving group class

Acceptance criteria:

- Every donor template is valid SMARTS.
- Audit reports per-template match count and map-retention on all test
  molecules.
- No donor template is active in `flavonoid_reaction_family_templates.hdf5`
  until it passes the audit and a human decision promotes it.

## Phase 4: Add A Cross-Flavonoid Target Panel

Goal: stop optimizing only for hesperidin.

Create:

- `config/flavonoid_target_panel.csv`
- update `scripts/run_ablation_experiments.py` or add
  `scripts/run_panel_ablation_experiments.py`
- update `scripts/compare_ablation_results.py` or add
  `scripts/compare_panel_ablation_results.py`
- `outputs/panel_ablation/`

Suggested initial panel:

```text
target_name,target_class,expected_sugar_family,smiles
hesperidin,flavanone_glycoside,neohesperidose,...
naringin,flavanone_glycoside,neohesperidose,...
neohesperidin,flavanone_glycoside,neohesperidose,...
narirutin,flavanone_glycoside,rutinose,...
rutin,flavonol_glycoside,rutinose,...
quercitrin,flavonol_glycoside,rhamnoside,...
```

If exact SMILES are introduced, validate them with RDKit and record source or
generation method. Do not silently use low-confidence target structures.

Acceptance criteria:

- The runner can execute the current nine ablations for one target or for the
  whole panel.
- The comparison report is target-aware and includes:
  - per-target strict/trusted/bridge/ZINC solved counts
  - aggregate solved counts
  - sugar-family breakdown
  - bridge-closed vs non-virtual separation
- Hesperidin-only output remains reproducible or is clearly superseded.

## Phase 5: Improve Reporting And Scoring Boundaries

Goal: make it impossible to confuse diagnostic bridge closure with real solved
routes.

Update reports to include these fields per route:

- `route_validity_class`:
  - `strict_trusted_solved`
  - `bridge_closed_connectivity`
  - `zinc_baseline_solved`
  - `unsolved`
- `uses_virtual_bridge`
- `uses_sugar_gap_bridge`
- `uses_aromatic_glycoside_leaf`
- `uses_donor_surrogate`
- `contains_protected_sugar_artifact`
- `map_retention_min`

Also fix or relabel the current `sugar donor` statistic in
`compare_ablation_results.py`. Right now it is inferred from glycoside/sugar
keywords and does not prove donor-layer chemistry.

Acceptance criteria:

- `outputs/ablation/ablation_report.md` and JSON distinguish route classes.
- The summary cannot imply that bridge closure proves strict buyability.
- Donor-surrogate routes, if any, are reported separately from neutral sugar
  bridge routes.

## Recommended Execution Order

1. Phase 1: evidence review scaffold and report.
2. Phase 2: linkage identity cleanup for rutinose/neohesperidose.
3. Re-run experiments 8 and 9 plus comparison to confirm bridge behavior.
4. Phase 5 reporting cleanup, especially route validity class.
5. Phase 3 donor-surrogate sandbox, disabled by default.
6. Phase 4 target panel after the donor/bridge semantics are stable.

## Commands For Regression Checks

Run these after each phase that changes scripts, templates, stock, or configs:

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_reaction_family_templates.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_ablation_experiments.py --experiments 8,9
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_ablation_results.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/generate_route_gap_worklist.py
/home/ljx/miniforge3/envs/retro/bin/python scripts/build_sugar_bridge_layer.py
```

Before claiming final improvement, run the full matrix:

```bash
/home/ljx/miniforge3/envs/retro/bin/python scripts/run_ablation_experiments.py --experiments 1,2,3,4,5,6,7,8,9
/home/ljx/miniforge3/envs/retro/bin/python scripts/compare_ablation_results.py
```

## Final Deliverables

At minimum, deliver:

- Updated/new scripts with conservative defaults.
- New evidence review CSV and audit markdown.
- Updated metadata only where evidence justifies it.
- Updated ablation report and structured JSON.
- A short final note explaining:
  - what changed
  - whether strict/trusted solved changed from zero
  - whether any bridge rows were rejected or promoted
  - whether donor templates remain disabled or were promoted
  - exact commands run and any failures

Your final reply must be in Chinese.