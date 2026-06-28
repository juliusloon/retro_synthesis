#!/usr/bin/env python3
"""
验证黄酮类模板在橙皮苷上的匹配和裂解行为
改进版：合并模板元数据以显示 reaction_id/name，并更鲁棒地提取裂解产物。
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

TEMPLATES_HDF5 = Path("/home/ljx/retro_synthesis/templates/flavonoid_templates.hdf5")
METADATA_CSV = Path("/home/ljx/retro_synthesis/templates/flavonoid_templates_metadata.csv")
STOCK_META = Path("/home/ljx/retro_synthesis/templates/flavonoid_stock_metadata.csv")
TEMPLATE_DIR = Path("/home/ljx/retro_synthesis/templates")

HESPERIDIN_SMILES = "COc1ccc(C2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)[C@H](O)[C@@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O"


def prepare_stock_if_needed():
    if STOCK_META.exists():
        return
    print("stock元数据不存在，正在生成...")
    files = {
        "flavonoid_precursors": TEMPLATE_DIR / "stock_flavonoid_precursors.csv",
        "glycosyl_donors": TEMPLATE_DIR / "stock_glycosyl_donors.csv",
        "protecting_reagents": TEMPLATE_DIR / "stock_protecting_reagents.csv",
        "catalysts_reagents": TEMPLATE_DIR / "stock_catalysts_reagents.csv",
    }
    metadata = []
    for category, fpath in files.items():
        try:
            df = pd.read_csv(fpath)
        except Exception as e:
            print(f"  警告：无法读取 {fpath.name}: {e}")
            continue
        smiles_col = 'smiles' if 'smiles' in df.columns else [c for c in df.columns if 'smiles' in c.lower()][0]
        for _, row in df.iterrows():
            smi = str(row[smiles_col]).strip()
            if smi and smi.lower() != 'nan':
                metadata.append({
                    'smiles': smi,
                    'category': category,
                    'name': row.get('name', ''),
                    'source_file': fpath.name,
                })
    meta_df = pd.DataFrame(metadata).drop_duplicates(subset=['smiles'])
    meta_df.to_csv(STOCK_META, index=False)
    print(f"已生成stock元数据: {STOCK_META} ({len(meta_df)} 条)")


def load_stock_smiles():
    df = pd.read_csv(STOCK_META)
    return set(df['smiles'].astype(str).str.strip())


def mol_to_smiles_safe(mol):
    """尽可能把反应产物分子转成 canonical SMILES，失败返回空字符串。"""
    try:
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        try:
            mol.UpdatePropertyCache(strict=False)
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        except Exception:
            return ""


def make_result(reaction_id, name, classification, priority, retro_template,
                matched, cracked, match_count, notes, n_product_sets=0,
                products_in_stock=False, product_smiles=""):
    return {
        'reaction_id': reaction_id,
        'name': name,
        'classification': classification,
        'priority': priority,
        'retro_template': retro_template,
        'matched': matched,
        'cracked': cracked,
        'match_count': match_count,
        'n_product_sets': n_product_sets,
        'notes': notes,
        'products_in_stock': products_in_stock,
        'product_smiles': product_smiles,
    }


def main():
    print("=" * 70)
    print("黄酮类模板对橙皮苷的验证")
    print("=" * 70)

    prepare_stock_if_needed()
    df = pd.read_hdf(TEMPLATES_HDF5, key='table')
    print(f"\n模板总数: {len(df)}")

    # 加载元数据并与 HDF5 按行索引对齐
    if METADATA_CSV.exists():
        meta = pd.read_csv(METADATA_CSV)
        meta = meta.set_index('template_code')
    else:
        meta = pd.DataFrame(index=range(len(df)))
    # 确保索引对齐
    meta = meta.reindex(range(len(df)))

    target = Chem.MolFromSmiles(HESPERIDIN_SMILES)
    if target is None:
        print("错误：无法解析橙皮苷SMILES")
        return
    print(f"橙皮苷原子数: {target.GetNumAtoms()}")
    print(f"橙皮苷重原子数: {target.GetNumHeavyAtoms()}")

    stock_set = load_stock_smiles()
    print(f"专用stock中分子数: {len(stock_set)}")

    results = []
    matched_count = 0
    cracked_count = 0

    for idx, row in df.iterrows():
        retro_template = row['retro_template']
        meta_row = meta.loc[idx] if idx in meta.index else {}
        reaction_id = meta_row.get('reaction_id', f'template_{idx}')
        name = meta_row.get('name', '')
        classification = meta_row.get('classification', row.get('classification', ''))
        priority = meta_row.get('priority', '')

        try:
            rxn = AllChem.ReactionFromSmarts(retro_template)
            if rxn is None:
                results.append(make_result(reaction_id, name, classification, priority, retro_template,
                                           False, False, 0, 'ReactionFromSmarts None'))
                continue

            n_react = rxn.GetNumReactantTemplates()
            if n_react == 0:
                results.append(make_result(reaction_id, name, classification, priority, retro_template,
                                           False, False, 0, 'no reactant template'))
                continue

            matches = target.GetSubstructMatches(rxn.GetReactantTemplate(0))
            match_count = len(matches)

            if match_count == 0:
                results.append(make_result(reaction_id, name, classification, priority, retro_template,
                                           False, False, 0, 'no substructure match'))
                continue

            matched_count += 1

            try:
                products_list = rxn.RunReactants((target,))
                n_product_sets = len(products_list)
            except Exception as e:
                results.append(make_result(reaction_id, name, classification, priority, retro_template,
                                           True, False, match_count, f'RunReactants: {str(e)[:80]}'))
                continue

            if n_product_sets == 0:
                results.append(make_result(reaction_id, name, classification, priority, retro_template,
                                           True, False, match_count, 'matched but no products'))
                continue

            cracked_count += 1

            product_smiles_list = []
            all_in_stock = True
            for products in products_list[:5]:
                for p in products:
                    p_smi = mol_to_smiles_safe(p)
                    if not p_smi:
                        all_in_stock = False
                        continue
                    p_mol = Chem.MolFromSmiles(p_smi)
                    canonical = Chem.MolToSmiles(p_mol, isomericSmiles=True) if p_mol else p_smi
                    product_smiles_list.append(canonical)
                    if canonical not in stock_set:
                        all_in_stock = False

            unique_products = list(dict.fromkeys(product_smiles_list))[:3]
            product_str = "; ".join(unique_products)
            notes = f'products: {product_str[:180]}'

            results.append(make_result(reaction_id, name, classification, priority, retro_template,
                                       True, True, match_count, notes, n_product_sets,
                                       all_in_stock, product_str))

        except Exception as e:
            results.append(make_result(reaction_id, name, classification, priority, retro_template,
                                       False, False, 0, f'exception: {str(e)[:80]}'))

    df_results = pd.DataFrame(results)

    print(f"\n匹配到橙皮苷的模板: {matched_count}")
    print(f"能裂解出前体的模板: {cracked_count}")
    print(f"裂解产物全在stock中: {df_results['products_in_stock'].sum()}")

    print(f"\n=== 按分类统计匹配情况 ===")
    summary = df_results.groupby('classification').agg({
        'matched': 'sum',
        'cracked': 'sum',
        'products_in_stock': 'sum',
        'reaction_id': 'count'
    }).rename(columns={'reaction_id': 'total'})
    print(summary.to_string())

    print(f"\n=== 能裂解的模板列表 ===")
    cracked_df = df_results[df_results['cracked']].copy()
    print(cracked_df[['reaction_id', 'name', 'classification', 'priority', 'match_count',
                      'n_product_sets', 'products_in_stock', 'product_smiles']].to_string())

    output_path = Path("/home/ljx/retro_synthesis/templates/template_validation_on_hesperidin.csv")
    df_results.to_csv(output_path, index=False)
    print(f"\n已保存验证报告: {output_path}")

    print(f"\n=== 关键橙皮苷相关模板状态 ===")
    key_templates = ['GLY_009', 'GLY_023', 'CYC_015', 'DEP_001', 'PROT_004', 'GLY_021']
    for rid in key_templates:
        subset = df_results[df_results['reaction_id'] == rid]
        if len(subset) > 0:
            r = subset.iloc[0]
            print(f"{rid}: matched={r['matched']}, cracked={r['cracked']}, stock={r['products_in_stock']}, notes={r['notes'][:120]}")
        else:
            print(f"{rid}: 不在成功模板中（可能在失败列表里）")

    # 额外输出失败的关键模板信息
    failed_path = Path("/home/ljx/retro_synthesis/templates/flavonoid_templates_failed.csv")
    if failed_path.exists():
        failed_df = pd.read_csv(failed_path)
        print(f"\n=== 失败列表中的关键模板 ===")
        for rid in key_templates:
            subset = failed_df[failed_df['reaction_id'] == rid]
            if len(subset) > 0:
                r = subset.iloc[0]
                print(f"{rid}: 失败原因={r['error']}")


if __name__ == '__main__':
    main()
