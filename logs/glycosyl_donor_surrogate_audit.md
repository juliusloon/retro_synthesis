# Glycosyl Donor Surrogate Audit

Generated: 2026-07-01 11:02

## Summary

- Total donor templates: 7
- Valid SMARTS: 7
- Atom mapped templates: 2
- Active expansion: 0

## Template Details

### glycosyl_bromide_donor

- **Leaving group**: bromide
- **Valid SMARTS**: True
- **Atom mapped**: False
- **Active expansion**: False
- **Template status**: placeholder_unmapped
- **Validation status**: failed_map_retention
- **Audit reason**: disabled until mapped donor identity and test molecules pass
- **Average map retention**: 0.00
- **Min map retention**: 0.00

**Test molecule results:**

| Molecule | Map Retention | Expected Donor Found | Donor InChIKeys |
|----------|---------------|---:|---|
| hesperidin | 0.00 |  | `` |
| naringin | 0.00 |  | `` |
| rutin | 0.00 |  | `` |
| phenyl_O_glucoside | 0.00 |  | `` |

### glycosyl_chloride_donor

- **Leaving group**: chloride
- **Valid SMARTS**: True
- **Atom mapped**: False
- **Active expansion**: False
- **Template status**: placeholder_unmapped
- **Validation status**: failed_map_retention
- **Audit reason**: disabled until mapped donor identity and test molecules pass
- **Average map retention**: 0.00
- **Min map retention**: 0.00

**Test molecule results:**

| Molecule | Map Retention | Expected Donor Found | Donor InChIKeys |
|----------|---------------|---:|---|
| hesperidin | 0.00 |  | `` |
| naringin | 0.00 |  | `` |
| rutin | 0.00 |  | `` |
| phenyl_O_glucoside | 0.00 |  | `` |

### beta_rutinosyl_chloride_cw_donor_sandbox

- **Leaving group**: chloride
- **Valid SMARTS**: True
- **Atom mapped**: True
- **Active expansion**: False
- **Template status**: atom_mapped_sandbox
- **Validation status**: canonical_hesperidin_exact_donor_pass_reverse_match_template_required_for_narirutin_rutin
- **Audit reason**: inactive sandbox for confirmed CW beta rutinosyl chloride donor; after canonicalizing input SMILES it produces the expected donor for hesperidin, while narirutin/rutin require the paired reverse-match sandbox because RDKit reaction chirality is match-direction dependent
- **Average map retention**: 0.67
- **Min map retention**: 0.00

- **Expected donor**: hexa-O-acetyl-beta-rutinosyl chloride
- **Expected donor InChIKey**: `BQIWKKDEVVDAAJ-FXXSWGFCSA-N`

**Test molecule results:**

| Molecule | Map Retention | Expected Donor Found | Donor InChIKeys |
|----------|---------------|---:|---|
| hesperidin_local | 1.00 | True | `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` |
| hesperidin_pubchem | 1.00 | True | `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` |
| narirutin_pubchem | 1.00 | False | `BQIWKKDEVVDAAJ-IRQRQJHYSA-N` |
| rutin_pubchem | 1.00 | False | `BQIWKKDEVVDAAJ-IRQRQJHYSA-N` |
| naringin_pubchem | 0.00 | False | `` |
| neohesperidin_pubchem | 0.00 | False | `` |

**Example precursor sets:**

- **hesperidin_local**: `COc1ccc(C2CC(=O)c3c(O)cc(O)cc3O2)cc1O.CC(=O)O[C@@H]1[C@@H](OC(C)=O)[C@H](Cl)O[C@H](CO[C@@H]2O[C@@H](C)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@H]2OC(C)=O)[C@H]1OC(C)=O`
- **hesperidin_pubchem**: `COc1ccc([C@@H]2CC(=O)c3c(O)cc(O)cc3O2)cc1O.CC(=O)O[C@@H]1[C@@H](OC(C)=O)[C@H](Cl)O[C@H](CO[C@@H]2O[C@@H](C)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@H]2OC(C)=O)[C@H]1OC(C)=O`
- **narirutin_pubchem**: `O=C1C[C@@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21.CC(=O)O[C@@H]1[C@@H](OC(C)=O)[C@@H](Cl)O[C@H](CO[C@@H]2O[C@@H](C)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@H]2OC(C)=O)[C@@H]1OC(C)=O`
- **rutin_pubchem**: `O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12.CC(=O)O[C@@H]1[C@@H](OC(C)=O)[C@@H](Cl)O[C@H](CO[C@@H]2O[C@@H](C)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@H]2OC(C)=O)[C@@H]1OC(C)=O`
- **naringin_pubchem**: `N/A`
- **neohesperidin_pubchem**: `N/A`

### beta_rutinosyl_chloride_cw_donor_reverse_match_sandbox

- **Leaving group**: chloride
- **Valid SMARTS**: True
- **Atom mapped**: True
- **Active expansion**: False
- **Template status**: atom_mapped_sandbox
- **Validation status**: narirutin_rutin_exact_donor_pass_hesperidin_covered_by_paired_template
- **Audit reason**: paired inactive sandbox for the opposite RDKit substructure-match direction; flips product-side map 3 and map 22 chirality so canonical narirutin/rutin produce the same confirmed CW donor InChIKey
- **Average map retention**: 0.67
- **Min map retention**: 0.00

- **Expected donor**: hexa-O-acetyl-beta-rutinosyl chloride
- **Expected donor InChIKey**: `BQIWKKDEVVDAAJ-FXXSWGFCSA-N`

**Test molecule results:**

| Molecule | Map Retention | Expected Donor Found | Donor InChIKeys |
|----------|---------------|---:|---|
| hesperidin_local | 1.00 | False | `BQIWKKDEVVDAAJ-IRQRQJHYSA-N` |
| hesperidin_pubchem | 1.00 | False | `BQIWKKDEVVDAAJ-IRQRQJHYSA-N` |
| narirutin_pubchem | 1.00 | True | `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` |
| rutin_pubchem | 1.00 | True | `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` |
| naringin_pubchem | 0.00 | False | `` |
| neohesperidin_pubchem | 0.00 | False | `` |

**Example precursor sets:**

- **hesperidin_local**: `COc1ccc(C2CC(=O)c3c(O)cc(O)cc3O2)cc1O.CC(=O)O[C@@H]1[C@@H](OC(C)=O)[C@@H](Cl)O[C@H](CO[C@@H]2O[C@@H](C)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@H]2OC(C)=O)[C@@H]1OC(C)=O`
- **hesperidin_pubchem**: `COc1ccc([C@@H]2CC(=O)c3c(O)cc(O)cc3O2)cc1O.CC(=O)O[C@@H]1[C@@H](OC(C)=O)[C@@H](Cl)O[C@H](CO[C@@H]2O[C@@H](C)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@H]2OC(C)=O)[C@@H]1OC(C)=O`
- **narirutin_pubchem**: `O=C1C[C@@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21.CC(=O)O[C@@H]1[C@@H](OC(C)=O)[C@H](Cl)O[C@H](CO[C@@H]2O[C@@H](C)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@H]2OC(C)=O)[C@H]1OC(C)=O`
- **rutin_pubchem**: `O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12.CC(=O)O[C@@H]1[C@@H](OC(C)=O)[C@H](Cl)O[C@H](CO[C@@H]2O[C@@H](C)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@H]2OC(C)=O)[C@H]1OC(C)=O`
- **naringin_pubchem**: `N/A`
- **neohesperidin_pubchem**: `N/A`

### glycosyl_trichloroacetimidate_donor

- **Leaving group**: trichloroacetimidate
- **Valid SMARTS**: True
- **Atom mapped**: False
- **Active expansion**: False
- **Template status**: placeholder_unmapped
- **Validation status**: failed_map_retention
- **Audit reason**: disabled until mapped donor identity and test molecules pass
- **Average map retention**: 0.00
- **Min map retention**: 0.00

**Test molecule results:**

| Molecule | Map Retention | Expected Donor Found | Donor InChIKeys |
|----------|---------------|---:|---|
| hesperidin | 0.00 |  | `` |
| naringin | 0.00 |  | `` |
| rutin | 0.00 |  | `` |
| phenyl_O_glucoside | 0.00 |  | `` |

### thioglycoside_donor

- **Leaving group**: thiomethyl
- **Valid SMARTS**: True
- **Atom mapped**: False
- **Active expansion**: False
- **Template status**: placeholder_unmapped
- **Validation status**: failed_map_retention
- **Audit reason**: disabled until mapped donor identity and test molecules pass
- **Average map retention**: 0.00
- **Min map retention**: 0.00

**Test molecule results:**

| Molecule | Map Retention | Expected Donor Found | Donor InChIKeys |
|----------|---------------|---:|---|
| hesperidin | 0.00 |  | `` |
| naringin | 0.00 |  | `` |
| rutin | 0.00 |  | `` |
| phenyl_O_glucoside | 0.00 |  | `` |

### glycosyl_acetate_surrogate

- **Leaving group**: acetate
- **Valid SMARTS**: True
- **Atom mapped**: False
- **Active expansion**: False
- **Template status**: placeholder_unmapped
- **Validation status**: failed_map_retention
- **Audit reason**: disabled until mapped donor identity and test molecules pass; surrogate only - clearly mark as surrogate
- **Average map retention**: 0.00
- **Min map retention**: 0.00

**Test molecule results:**

| Molecule | Map Retention | Expected Donor Found | Donor InChIKeys |
|----------|---------------|---:|---|
| hesperidin | 0.00 |  | `` |
| naringin | 0.00 |  | `` |
| rutin | 0.00 |  | `` |
| phenyl_O_glucoside | 0.00 |  | `` |

## Expected Donor Coverage Across Templates

| Expected Donor InChIKey | Molecule | Found By Templates |
|---|---|---|
| `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` | hesperidin_local | beta_rutinosyl_chloride_cw_donor_sandbox |
| `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` | hesperidin_pubchem | beta_rutinosyl_chloride_cw_donor_sandbox |
| `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` | naringin_pubchem | - |
| `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` | narirutin_pubchem | beta_rutinosyl_chloride_cw_donor_reverse_match_sandbox |
| `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` | neohesperidin_pubchem | - |
| `BQIWKKDEVVDAAJ-FXXSWGFCSA-N` | rutin_pubchem | beta_rutinosyl_chloride_cw_donor_reverse_match_sandbox |

This table is the panel-level interpretation for paired sandbox templates. A single template is match-direction specific; coverage can be achieved by the paired inactive sandbox templates without promoting either to production.

## Acceptance Criteria

- [x] Every donor template is valid SMARTS: **YES**
- [ ] Every donor template is atom mapped: **NO**
- [x] Audit reports per-template map-retention: **YES**
- [x] Audit checks expected donor identity when configured: **YES**
- [x] No donor template is active: **YES**
- [ ] All templates pass min_retained_map_ratio >= 0.8: **NO**

## Notes

- All templates are disabled by default (active_expansion=false)
- Test molecules are RDKit-canonicalized before reaction application to remove raw SMILES atom-order artifacts
- Current donor templates are placeholders until atom-mapped donor identity is encoded
- The CW beta rutinosyl chloride entry is an inactive sandbox, not a production expansion template
- Templates must pass audit and human decision before promotion
- Do not enable in flavonoid_reaction_family_templates.hdf5 until validated
