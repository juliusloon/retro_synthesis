# Glycosyl Donor Surrogate Audit

Generated: 2026-06-30 19:59

## Summary

- Total donor templates: 5
- Valid SMARTS: 5
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

| Molecule | Map Retention |
|----------|---------------|
| hesperidin | 0.00 |
| naringin | 0.00 |
| rutin | 0.00 |
| phenyl_O_glucoside | 0.00 |

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

| Molecule | Map Retention |
|----------|---------------|
| hesperidin | 0.00 |
| naringin | 0.00 |
| rutin | 0.00 |
| phenyl_O_glucoside | 0.00 |

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

| Molecule | Map Retention |
|----------|---------------|
| hesperidin | 0.00 |
| naringin | 0.00 |
| rutin | 0.00 |
| phenyl_O_glucoside | 0.00 |

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

| Molecule | Map Retention |
|----------|---------------|
| hesperidin | 0.00 |
| naringin | 0.00 |
| rutin | 0.00 |
| phenyl_O_glucoside | 0.00 |

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

| Molecule | Map Retention |
|----------|---------------|
| hesperidin | 0.00 |
| naringin | 0.00 |
| rutin | 0.00 |
| phenyl_O_glucoside | 0.00 |

## Acceptance Criteria

- [x] Every donor template is valid SMARTS: **YES**
- [ ] Every donor template is atom mapped: **NO**
- [x] Audit reports per-template map-retention: **YES**
- [x] No donor template is active: **YES**
- [ ] All templates pass min_retained_map_ratio >= 0.8: **NO**

## Notes

- All templates are disabled by default (active_expansion=false)
- Current donor templates are placeholders until atom-mapped donor identity is encoded
- Templates must pass audit and human decision before promotion
- Do not enable in flavonoid_reaction_family_templates.hdf5 until validated
