# Flavonoid Literature Template Build Report

- Generated at: 2026-06-28T12:12:44.167590+00:00
- Source directory: `/home/ljx/retro_synthesis/files`
- Extracted reaction entries: 508
- RDKit-valid retro templates: 508
- Unique full templates: 508
- Conservative AiZynthFinder-ready templates: 508
- Entries marked AiZynthFinder-ready: 508
- Entries with nearby condition text: 437
- Entries containing generic/wildcard groups: 0

## Classification Counts

- chalcone_synthesis: 115
- literature_reaction: 80
- biosynthesis: 78
- flavonoid_synthesis: 51
- isoflavone_synthesis: 45
- protection: 40
- flavanone_synthesis: 31
- reduction: 17
- oxidation: 15
- deprotection: 9
- flavonoid_cyclization: 6
- glycosylation: 6
- methylation: 6
- cyclization: 5
- c_glycosylation: 2
- o_glycosylation: 2

## Generic Resolution Counts

- none: 0

## Source Counts

- 26-06-26 黄酮化合物的合成.md: 404
- 26-06-18 黄酮化合物的生物合成.md: 50
- 26-06-25 黄酮化合物的性质.md: 39
- 26-06-21 黄酮类化合物立体化学.md: 8
- 26-06-15 黄酮化合物的结构.md: 4
- 26-06-21 黄酮化合物研究进展（1）：黄酮类化合物.md: 3

## Outputs

- all_reactants_inchikeys: `/home/ljx/retro_synthesis/templates/literature_curated/flavonoid_literature_all_reactants_stock_inchikeys.txt`
- bio_hdf5: `/home/ljx/retro_synthesis/templates/flavonoid_biosynthesis/flavonoid_biosynthesis_templates.hdf5`
- candidates_csv: `/home/ljx/retro_synthesis/templates/literature_curated/flavonoid_literature_stock_candidates.csv`
- csv_gz: `/home/ljx/retro_synthesis/templates/literature_curated/flavonoid_literature_templates.csv.gz`
- curated_json: `/home/ljx/retro_synthesis/templates/literature_curated/flavonoid_literature_reactions.json`
- generic_review_csv: `/home/ljx/retro_synthesis/templates/literature_curated/generic_group_review.csv`
- generic_review_json: `/home/ljx/retro_synthesis/templates/literature_curated/generic_group_review.json`
- generic_review_md: `/home/ljx/retro_synthesis/templates/literature_curated/generic_group_review.md`
- hdf5: `/home/ljx/retro_synthesis/templates/literature_curated/flavonoid_literature_templates.hdf5`
- invalid_csv: `/home/ljx/retro_synthesis/templates/literature_curated/flavonoid_literature_invalid_reactions.csv`
- manual_review_override_summary: `/home/ljx/retro_synthesis/templates/literature_curated/manual_review_override_summary.json`
- metadata_csv: `/home/ljx/retro_synthesis/templates/literature_curated/flavonoid_literature_templates_metadata.csv`
- ready_csv_gz: `/home/ljx/retro_synthesis/templates/literature_curated/flavonoid_literature_templates_aizynthfinder_ready.csv.gz`
- ready_hdf5: `/home/ljx/retro_synthesis/templates/literature_curated/flavonoid_literature_templates_aizynthfinder_ready.hdf5`
- root_csv_gz: `/home/ljx/retro_synthesis/templates/flavonoid_templates.csv.gz`
- root_hdf5: `/home/ljx/retro_synthesis/templates/flavonoid_templates.hdf5`
- root_stock: `/home/ljx/retro_synthesis/templates/flavonoid_stock_inchikeys.txt`
- strict_inchikeys: `/home/ljx/retro_synthesis/templates/literature_curated/flavonoid_literature_stock_inchikeys.txt`
- structural_hdf5: `/home/ljx/retro_synthesis/templates/flavonoid_structural_diversity/flavonoid_structural_templates.hdf5`

## Curation Notes

- Manual review overrides from `generic_group_review.md` are applied before export: finite placeholders are expanded to concrete structures, enzyme/cofactor placeholders are moved to conditions, and broad aryl/alkyl placeholders are collapsed to parent skeleton templates.
- The default `/home/ljx/retro_synthesis/templates/flavonoid_templates.hdf5` now includes every RDKit-valid curated entry; no generic placeholders remain in the ready export.
- Full reviewed entries, including placeholder annotations, manual review notes, and possible structures, remain in `literature_curated/flavonoid_literature_reactions.json`; `generic_group_review.*` is regenerated as an empty residual review set when no generic entries remain.
- `literature_curated/manual_review_override_summary.json` records how many variants each manually reviewed parent reaction generated.
- One truncated source SMILES (`[I]c1ccc(OC(=O)C)`) is repaired to `[I]c1ccc(OC(=O)C)cc1` during extraction so the corresponding chalcone step remains usable.
- The strict stock file includes valid reactant molecules that do not also appear as products in the extracted reaction set, to avoid marking route intermediates as purchasable by default.
- `flavonoid_literature_reactions.json` keeps source context and extracted condition text for manual audit.
