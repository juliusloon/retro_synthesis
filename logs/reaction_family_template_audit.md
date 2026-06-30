# Reaction Family Template Audit

Generated: 2026-06-30 19:28

- Source: `/home/ljx/retro_synthesis/config/reaction_families.json`
- Output HDF5: `/home/ljx/retro_synthesis/templates/reaction_families/flavonoid_reaction_family_templates.hdf5`
- Output CSV: `/home/ljx/retro_synthesis/templates/reaction_families/flavonoid_reaction_family_templates.csv`

- Audit CSV: `/home/ljx/retro_synthesis/templates/reaction_families/flavonoid_reaction_family_templates_audit.csv`
- Families in source: 11
- Active templates written to search library: 3
- Inactive templates retained only for audit: 8

## Test Molecules

- **hesperidin**: `COc1ccc(C2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)[C@H](O)[C@@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O`
- **hesperetin**: `COc1ccc(C2CC(=O)c3c(O)cc(O)cc3O2)cc1O`
- **naringenin**: `OC1CC(=O)c2c(O)cc(O)cc2O1`
- **hydroxychalcone_2**: `O=C(c1ccccc1)/C=C/c1ccccc1O`
- **simple_chalcone**: `O=C(c1ccccc1)/C=C/c1ccccc1`

## Template Summary

| Family | Active | Classification | Template Class | Valid SMARTS | Max Map Retention | Passes Min Ratio | Evidence Family |
|---|---:|---|---|---:|---:|---:|---|
| o_glycoside_cleavage_pyranose | ✓ | o_glycoside_cleavage | mechanistic_family_template | True | 0.814 | ✓ | aryl_O_glycoside_disconnection |
| aryl_methyl_ether_cleavage | ✓ | o_demethylation | mechanistic_family_template | True | 0.930 | ✓ | aryl_OMe_demethylation |
| hydroxy_flavanone_to_chalcone | ✗ | flavanone_chalcone_opening | mechanistic_family_template | True | 0.558 | ✗ | flavanone_to_chalcone |
| chalcone_retro_aldol | ✗ | chalcone_retro_aldol | mechanistic_family_template | True | 0.000 | ✗ | chalcone_retro_aldol |
| activated_glycosyl_donor_disconnection | ✗ | activated_glycosyl_donor | mechanistic_family_template | True | 0.814 | ✓ | activated_glycosyl_donor_disconnection |
| phenol_acetylation | ✗ | phenol_protection | mechanistic_family_template | True | 0.000 | ✗ | phenol_protection |
| phenol_deacetylation | ✗ | phenol_deprotection | mechanistic_family_template | True | 0.000 | ✗ | phenol_deprotection |
| phenol_benzylation | ✗ | phenol_protection | mechanistic_family_template | True | 0.000 | ✗ | phenol_protection |
| phenol_debenzylation | ✗ | phenol_deprotection | mechanistic_family_template | True | 0.000 | ✗ | phenol_deprotection |
| aryl_O_methylation | ✓ | o_methylation | mechanistic_family_template | True | 0.930 | ✓ | aryl_O_methylation |
| chalcone_aldol_condensation | ✗ | chalcone_aldol | mechanistic_family_template | True | 0.000 | ✗ | chalcone_aldol_condensation |

## Detailed Validation Results

### o_glycoside_cleavage_pyranose

- **Classification**: o_glycoside_cleavage
- **Template Class**: mechanistic_family_template
- **Evidence Family**: aryl_O_glycoside_disconnection
- **Active Expansion**: True
- **Audit Reason**: 
- **Min Retained Map Ratio**: 0.8
- **Notes**: Disconnects aryl O-glycoside bonds to aglycone phenol plus sugar fragment.

| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |
|---|---|---|---:|---:|---|
| hesperidin | True | True | 1 | 0.814 | `Oc1[cH:14][c:12]([OH:13])[c:11]2[c:39]([cH:38]1)[O` |
| hesperetin | True | False | 0 | 0.000 | `N/A` |
| naringenin | True | False | 0 | 0.000 | `N/A` |
| hydroxychalcone_2 | True | False | 0 | 0.000 | `N/A` |
| simple_chalcone | True | False | 0 | 0.000 | `N/A` |

### aryl_methyl_ether_cleavage

- **Classification**: o_demethylation
- **Template Class**: mechanistic_family_template
- **Evidence Family**: aryl_OMe_demethylation
- **Active Expansion**: True
- **Audit Reason**: 
- **Min Retained Map Ratio**: 0.8
- **Notes**: Converts aryl methyl ether to phenol plus methyl iodide surrogate.

| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |
|---|---|---|---:|---:|---|
| hesperidin | True | True | 1 | 0.930 | `Oc1[cH:4][cH:5][c:6]([CH:7]2[CH2:8][C:9](=[O:10])[` |
| hesperetin | True | True | 1 | 0.864 | `Oc1[cH:4][cH:5][c:6]([CH:7]2[CH2:8][C:9](=[O:10])[` |
| naringenin | True | False | 0 | 0.000 | `N/A` |
| hydroxychalcone_2 | True | False | 0 | 0.000 | `N/A` |
| simple_chalcone | True | False | 0 | 0.000 | `N/A` |

### hydroxy_flavanone_to_chalcone

- **Classification**: flavanone_chalcone_opening
- **Template Class**: mechanistic_family_template
- **Evidence Family**: flavanone_to_chalcone
- **Active Expansion**: False
- **Audit Reason**: max_map_retention=0.558, below threshold 0.8
- **Min Retained Map Ratio**: 0.8
- **Notes**: Mapped version of flavanone ring opening to chalcone. Long substituent retention is audited post-application.

| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |
|---|---|---|---:|---:|---|
| hesperidin | True | True | 2 | 0.558 | `O=C(/C=C/c1ccc(O[CH3:1])c([OH:43])c1)c1cc(O[C@@H:1` |
| hesperetin | True | True | 2 | 0.136 | `O=C(/C=C/c1ccc(O[CH3:1])c([OH:22])c1)c1cc(O)cc([OH` |
| naringenin | True | False | 0 | 0.000 | `N/A` |
| hydroxychalcone_2 | True | False | 0 | 0.000 | `N/A` |
| simple_chalcone | True | False | 0 | 0.000 | `N/A` |

### chalcone_retro_aldol

- **Classification**: chalcone_retro_aldol
- **Template Class**: mechanistic_family_template
- **Evidence Family**: chalcone_retro_aldol
- **Active Expansion**: False
- **Audit Reason**: does not match chalcone test molecules with acceptable map retention; requires separate SMARTS rewrite before expansion use
- **Min Retained Map Ratio**: 0.8
- **Notes**: Retro-aldol disconnection of chalcone at alpha-beta bond to acetophenone + benzaldehyde.

| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |
|---|---|---|---:|---:|---|
| hesperidin | True | False | 0 | 0.000 | `N/A` |
| hesperetin | True | False | 0 | 0.000 | `N/A` |
| naringenin | True | False | 0 | 0.000 | `N/A` |
| hydroxychalcone_2 | True | False | 0 | 0.000 | `N/A` |
| simple_chalcone | True | False | 0 | 0.000 | `N/A` |

### activated_glycosyl_donor_disconnection

- **Classification**: activated_glycosyl_donor
- **Template Class**: mechanistic_family_template
- **Evidence Family**: activated_glycosyl_donor_disconnection
- **Active Expansion**: False
- **Audit Reason**: identical SMARTS to o_glycoside_cleavage_pyranose; donor-side leaving group not encoded. Treat as abstract donor, not for solved evaluation.
- **Min Retained Map Ratio**: 0.8
- **Notes**: Abstract donor-side glycoside disconnect. Same SMARTS as o_glycoside_cleavage but semantically for activated donors. Does NOT count toward solved.

| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |
|---|---|---|---:|---:|---|
| hesperidin | True | True | 1 | 0.814 | `Oc1[cH:14][c:12]([OH:13])[c:11]2[c:39]([cH:38]1)[O` |
| hesperetin | True | False | 0 | 0.000 | `N/A` |
| naringenin | True | False | 0 | 0.000 | `N/A` |
| hydroxychalcone_2 | True | False | 0 | 0.000 | `N/A` |
| simple_chalcone | True | False | 0 | 0.000 | `N/A` |

### phenol_acetylation

- **Classification**: phenol_protection
- **Template Class**: mechanistic_family_template
- **Evidence Family**: phenol_protection
- **Active Expansion**: False
- **Audit Reason**: max_map_retention=0.0, does not match any test molecule
- **Min Retained Map Ratio**: 0.8
- **Notes**: Retro-acetylation of phenol acetate. Audit only - no test molecule match.

| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |
|---|---|---|---:|---:|---|
| hesperidin | True | False | 0 | 0.000 | `N/A` |
| hesperetin | True | False | 0 | 0.000 | `N/A` |
| naringenin | True | False | 0 | 0.000 | `N/A` |
| hydroxychalcone_2 | True | False | 0 | 0.000 | `N/A` |
| simple_chalcone | True | False | 0 | 0.000 | `N/A` |

### phenol_deacetylation

- **Classification**: phenol_deprotection
- **Template Class**: mechanistic_family_template
- **Evidence Family**: phenol_deprotection
- **Active Expansion**: False
- **Audit Reason**: max_map_retention=0.0, does not match any test molecule
- **Min Retained Map Ratio**: 0.8
- **Notes**: Deacetylation of phenol acetate. Audit only - no test molecule match.

| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |
|---|---|---|---:|---:|---|
| hesperidin | True | False | 0 | 0.000 | `N/A` |
| hesperetin | True | False | 0 | 0.000 | `N/A` |
| naringenin | True | False | 0 | 0.000 | `N/A` |
| hydroxychalcone_2 | True | False | 0 | 0.000 | `N/A` |
| simple_chalcone | True | False | 0 | 0.000 | `N/A` |

### phenol_benzylation

- **Classification**: phenol_protection
- **Template Class**: mechanistic_family_template
- **Evidence Family**: phenol_protection
- **Active Expansion**: False
- **Audit Reason**: max_map_retention=0.0, does not match any test molecule
- **Min Retained Map Ratio**: 0.8
- **Notes**: Retro-benzylation of phenol benzyl ether. Audit only - no test molecule match.

| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |
|---|---|---|---:|---:|---|
| hesperidin | True | False | 0 | 0.000 | `N/A` |
| hesperetin | True | False | 0 | 0.000 | `N/A` |
| naringenin | True | False | 0 | 0.000 | `N/A` |
| hydroxychalcone_2 | True | False | 0 | 0.000 | `N/A` |
| simple_chalcone | True | False | 0 | 0.000 | `N/A` |

### phenol_debenzylation

- **Classification**: phenol_deprotection
- **Template Class**: mechanistic_family_template
- **Evidence Family**: phenol_deprotection
- **Active Expansion**: False
- **Audit Reason**: max_map_retention=0.0, does not match any test molecule
- **Min Retained Map Ratio**: 0.8
- **Notes**: Debenzylation of phenol benzyl ether. Audit only - no test molecule match.

| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |
|---|---|---|---:|---:|---|
| hesperidin | True | False | 0 | 0.000 | `N/A` |
| hesperetin | True | False | 0 | 0.000 | `N/A` |
| naringenin | True | False | 0 | 0.000 | `N/A` |
| hydroxychalcone_2 | True | False | 0 | 0.000 | `N/A` |
| simple_chalcone | True | False | 0 | 0.000 | `N/A` |

### aryl_O_methylation

- **Classification**: o_methylation
- **Template Class**: mechanistic_family_template
- **Evidence Family**: aryl_O_methylation
- **Active Expansion**: True
- **Audit Reason**: 
- **Min Retained Map Ratio**: 0.8
- **Notes**: Retro-methylation of aryl methyl ether to phenol plus methyl iodide.

| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |
|---|---|---|---:|---:|---|
| hesperidin | True | True | 1 | 0.930 | `Oc1[cH:4][cH:5][c:6]([CH:7]2[CH2:8][C:9](=[O:10])[` |
| hesperetin | True | True | 1 | 0.864 | `Oc1[cH:4][cH:5][c:6]([CH:7]2[CH2:8][C:9](=[O:10])[` |
| naringenin | True | False | 0 | 0.000 | `N/A` |
| hydroxychalcone_2 | True | False | 0 | 0.000 | `N/A` |
| simple_chalcone | True | False | 0 | 0.000 | `N/A` |

### chalcone_aldol_condensation

- **Classification**: chalcone_aldol
- **Template Class**: mechanistic_family_template
- **Evidence Family**: chalcone_aldol_condensation
- **Active Expansion**: False
- **Audit Reason**: same SMARTS as chalcone_retro_aldol; disabled until chalcone SMARTS preserves atom maps under retention audit
- **Min Retained Map Ratio**: 0.8
- **Notes**: Same retro as chalcone_retro_aldol. Aldol condensation reverse: acetophenone + benzaldehyde -> chalcone.

| Molecule | Valid SMARTS | Matches | Outcomes | Max Map Retention | Example Precursors |
|---|---|---|---:|---:|---|
| hesperidin | True | False | 0 | 0.000 | `N/A` |
| hesperetin | True | False | 0 | 0.000 | `N/A` |
| naringenin | True | False | 0 | 0.000 | `N/A` |
| hydroxychalcone_2 | True | False | 0 | 0.000 | `N/A` |
| simple_chalcone | True | False | 0 | 0.000 | `N/A` |

## Notes

Only templates with `active_expansion: true` are written to the HDF5 and CSV search libraries. Inactive templates are retained in the audit report and audit CSV for traceability.

The expansion strategy also gates applied actions with `min_retained_map_ratio`, so target-validation retention and post-application route retention are both visible.

All templates must have valid SMARTS and atom maps. Templates with 0 map retention are flagged for review.
