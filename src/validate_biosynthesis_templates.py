#!/usr/bin/env python3
"""
Validate flavonoid biosynthesis templates against representative targets.
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

TEMPLATES_HDF5 = Path("/home/ljx/retro_synthesis/templates/flavonoid_biosynthesis/flavonoid_biosynthesis_templates.hdf5")
METADATA_CSV = Path("/home/ljx/retro_synthesis/templates/flavonoid_biosynthesis/flavonoid_biosynthesis_templates_metadata.csv")
NEW_STOCK_CSV = Path("/home/ljx/retro_synthesis/templates/flavonoid_biosynthesis/flavonoid_biosynthesis_stock.csv")
EXISTING_STOCK_CSV = Path("/home/ljx/retro_synthesis/templates/flavonoid_stock_metadata.csv")
STRUCTURAL_STOCK_CSV = Path("/home/ljx/retro_synthesis/templates/flavonoid_structural_diversity/flavonoid_structural_stock.csv")

TARGETS = {
    "hesperidin": "COc1ccc(C2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)[C@H](O)[C@@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O",
    "rutin": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O[C@@H]3O[C@H](CO[C@@H]4O[C@@H](C)[C@H](O)[C@@H](O)[C@H]4O)[C@@H](O)[C@H](O)[C@H]3O)ccc12",
    "naringenin": "O=C1C[C@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21",
    "eriodictyol": "O=C1C[C@H](c2ccc(O)c(O)c2)Oc2cc(O)cc(O)c21",
    "dihydrokaempferol": "O=C1c2c(O)cc(O)cc2O[C@H](c2ccc(O)cc2)[C@H]1O",
    "taxifolin": "O=C1c2c(O)cc(O)cc2O[C@H](c2ccc(O)c(O)c2)[C@H]1O",
    "kaempferol": "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",
    "quercetin": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
    "apigenin": "O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",
    "luteolin": "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
    "genistein": "Oc1cc2c(=O)c(-c3ccc(O)cc3)coc2cc1O",
    "daidzein": "O=c1c(-c2ccc(O)cc2)coc2ccccc12",
    "formononetin": "COc1ccc(-c2coc3ccccc3c2=O)cc1",
    "medicarpin": "COc1ccc2c(c1)O[C@H]1c3ccc(O)cc3OC[C@@H]21",
    "rotenone": "C=C(C)[C@H]1Oc2ccc3C(=O)[C@H]4c5cc(OC)c(OC)cc5OC[C@H]4Oc3c2C1",
    "naringin_dihydrochalcone": "Oc1cc(O)cc(O)c1C(=O)CCc2ccc(O)cc2",
    "silybin": "COc1cc([C@@H]2Oc3c(O)cc([C@H]4Oc5cc(O)cc(O)c5C(=O)[C@@H]4O)cc3[C@@H]2CO)ccc1O",
    "resveratrol": "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1",
    "parent_chalcone": "O=C(/C=C/c1ccccc1)c1ccccc1",
    "parent_flavanone": "O=C1CC(c2ccccc2)Oc2ccccc21",
    "parent_flavone": "O=c1cc(-c2ccccc2)oc2ccccc12",
}


def load_stock_smiles():
    smiles = set()
    for p in [NEW_STOCK_CSV, EXISTING_STOCK_CSV, STRUCTURAL_STOCK_CSV]:
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p)
            col = 'smiles' if 'smiles' in df.columns else [c for c in df.columns if 'smiles' in c.lower()][0]
            smiles.update(df[col].astype(str).str.strip())
        except Exception as e:
            print(f"Warning reading {p}: {e}")
    return smiles


def mol_to_smiles_safe(mol):
    try:
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        try:
            mol.UpdatePropertyCache(strict=False)
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        except Exception:
            return ""


def validate_target(target_name, target_smiles, df_templates, meta, stock_set):
    target = Chem.MolFromSmiles(target_smiles)
    if target is None:
        print(f"  {target_name}: INVALID SMILES")
        return []
    results = []
    for idx, row in df_templates.iterrows():
        retro_template = row['retro_template']
        meta_row = meta.loc[idx] if idx in meta.index else {}
        reaction_id = meta_row.get('reaction_id', f'template_{idx}')
        name = meta_row.get('name', '')
        classification = meta_row.get('classification', row.get('classification', ''))
        priority = meta_row.get('priority', '')
        try:
            rxn = AllChem.ReactionFromSmarts(retro_template)
            if rxn is None or rxn.GetNumReactantTemplates() == 0:
                continue
            matches = target.GetSubstructMatches(rxn.GetReactantTemplate(0))
            if not matches:
                continue
            products_list = rxn.RunReactants((target,))
            n_product_sets = len(products_list)
            if n_product_sets == 0:
                continue
            product_smiles_list = []
            all_in_stock = True
            for products in products_list[:3]:
                for p in products:
                    p_smi = mol_to_smiles_safe(p)
                    if not p_smi:
                        all_in_stock = False
                        continue
                    product_smiles_list.append(p_smi)
                    if p_smi not in stock_set:
                        all_in_stock = False
            unique_products = list(dict.fromkeys(product_smiles_list))[:3]
            results.append({
                'target': target_name,
                'reaction_id': reaction_id,
                'name': name,
                'classification': classification,
                'priority': priority,
                'matched': True,
                'cracked': True,
                'match_count': len(matches),
                'n_product_sets': n_product_sets,
                'products_in_stock': all_in_stock,
                'product_smiles': "; ".join(unique_products),
            })
        except Exception as e:
            pass
    return results


def main():
    print("=" * 70)
    print("黄酮生物合成模板验证")
    print("=" * 70)

    df = pd.read_hdf(TEMPLATES_HDF5, key='table')
    print(f"模板总数: {len(df)}")

    if METADATA_CSV.exists():
        meta = pd.read_csv(METADATA_CSV)
        meta = meta.set_index('template_code')
    else:
        meta = pd.DataFrame(index=range(len(df)))
    meta = meta.reindex(range(len(df)))

    stock_set = load_stock_smiles()
    print(f"合并 stock 中分子数: {len(stock_set)}")

    all_results = []
    summary = []
    for name, smi in TARGETS.items():
        res = validate_target(name, smi, df, meta, stock_set)
        all_results.extend(res)
        summary.append({
            'target': name,
            'matched_templates': len(res),
            'cracked_templates': len(res),
            'products_in_stock_for_all': all(r['products_in_stock'] for r in res) if res else False,
        })
        print(f"  {name}: {len(res)} templates matched+&cracked")

    df_results = pd.DataFrame(all_results)
    out_csv = TEMPLATES_HDF5.parent / 'template_validation_biosynthesis.csv'
    df_results.to_csv(out_csv, index=False)
    print(f"\n已保存验证结果: {out_csv}")

    print("\n按目标汇总:")
    print(pd.DataFrame(summary).to_string(index=False))

    print("\n按分类匹配数:")
    if not df_results.empty:
        print(df_results.groupby('classification').size().to_string())


if __name__ == '__main__':
    main()
