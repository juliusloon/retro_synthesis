#!/usr/bin/env python3
"""
清洗并转换黄酮类专用模板为aizynthfinder格式
"""
import json
import re
import hashlib
import pandas as pd
from rdkit.Chem import AllChem
from pathlib import Path

TEMPLATE_FILES = [
    "templates_glycosylation.json",
    "templates_deprotection.json",
    "templates_cyclization.json",
    "templates_protection_synthesis.json",
    "templates_additional.json",
]

OUTPUT_DIR = Path("/home/ljx/retro_synthesis/templates")


def clean_smarts(smarts):
    """清洗SMARTS，修复常见非法写法"""
    if not smarts:
        return smarts
    
    original = smarts
    
    # 1. 修复苯环 [Ring1][=Branch1] 凯库勒式为芳香性表示
    # 模式A: 六个交替双键碳，无支链
    smarts = re.sub(
        r'\[C:(\d+)\]=\[C:(\d+)\]\[C:(\d+)\]=\[C:(\d+)\]\[C:(\d+)\]=\[C:(\d+)\]\[Ring1\]\[=Branch1\]',
        r'[c:\1]1[c:\2][c:\3][c:\4][c:\5][c:\6]1',
        smarts
    )
    
    # 模式B: 带前置支链的苯环，如 X-[C:a]=[C:b]...[Ring1][=Branch1]
    # 尝试识别 [C:a]=[C:b][C:c]=[C:d][C:e]=[C:f][Ring1][=Branch1] 并替换
    # 这里用更宽松的模式重复处理
    smarts = re.sub(
        r'\[C:(\d+)\]=\[C:(\d+)\]\[C:(\d+)\]=\[C:(\d+)\]\[C:(\d+)\]=\[C:(\d+)\]\[Ring1\]\[=Branch1\]',
        r'[c:\1]1[c:\2][c:\3][c:\4][c:\5][c:\6]1',
        smarts
    )
    
    # 2. 移除独立显式氢 [H:x]，但保留 [C@H:x] 等
    smarts = re.sub(r'(?<![A-Za-z@])\[H:(\d+)\]', '', smarts)
    
    # 3. 清理因移除氢导致的空括号、连续点号
    smarts = re.sub(r'\(\s*\)', '', smarts)
    smarts = re.sub(r'\.\.', '.', smarts)
    smarts = re.sub(r'\.\]', ']', smarts)
    smarts = re.sub(r'\]\(', '](', smarts)
    
    # 4. 修复 [Ring1][=Branch1] 残留（非苯环情况）
    # 简单替换为芳香小写c的环闭合可能不适用，先记录
    if '[Ring1][=Branch1]' in smarts:
        # 尝试通用替换：找到前面最近的 [C:n] 加环闭合
        pass  # 复杂情况后续处理
    
    return smarts


def reverse_reaction(smarts):
    """翻转反应方向：反应物>>产物 -> 产物>>反应物"""
    if '>>' not in smarts:
        return smarts
    left, right = smarts.split('>>', 1)
    return f"{right.strip()}>>{left.strip()}"


def validate_reaction(smarts):
    """验证反应SMARTS是否合法"""
    if '>>' not in smarts:
        return False, "no reaction arrow"
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        if rxn is None:
            return False, "None"
        n_react = rxn.GetNumReactantTemplates()
        n_prod = rxn.GetNumProductTemplates()
        if n_react == 0 or n_prod == 0:
            return False, f"reactants={n_react}, products={n_prod}"
        return True, f"reactants={n_react}, products={n_prod}"
    except Exception as e:
        return False, str(e)[:150]


def load_templates():
    templates_dir = OUTPUT_DIR
    all_templates = []
    for fname in TEMPLATE_FILES:
        fpath = templates_dir / fname
        with open(fpath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        items = data if isinstance(data, list) else data.get("templates", [])
        for item in items:
            item['source_file'] = fname
        all_templates.extend(items)
    return all_templates


def main():
    print("=" * 70)
    print("黄酮类模板清洗与转换")
    print("=" * 70)
    
    all_templates = load_templates()
    print(f"总计加载: {len(all_templates)} 个模板\n")
    
    rows = []
    failed = []
    
    for t in all_templates:
        rid = t.get('reaction_id', 'unknown')
        name = t.get('name', '')
        original = t.get('smarts', '')
        
        # 清洗
        cleaned = clean_smarts(original)
        
        # 翻转方向
        reversed_smarts = reverse_reaction(cleaned)
        
        # 验证
        ok, msg = validate_reaction(reversed_smarts)
        
        if ok:
            template_hash = hashlib.sha256(reversed_smarts.encode()).hexdigest()
            rows.append({
                'template_code': len(rows),
                'retro_template': reversed_smarts,
                'template_hash': template_hash,
                'classification': t.get('category', 'flavonoid_custom'),
                'library_occurence': 1,
                'reaction_id': rid,
                'name': name,
                'source_file': t.get('source_file', ''),
                'priority': t.get('priority', 'medium'),
                'doi': t.get('doi', ''),
                'conditions': t.get('conditions', ''),
                'original_smarts': original,
                'cleaned_smarts': cleaned,
            })
        else:
            failed.append({
                'reaction_id': rid,
                'name': name,
                'error': msg,
                'original': original,
                'cleaned': cleaned,
                'reversed': reversed_smarts,
            })
    
    print(f"成功: {len(rows)}")
    print(f"失败: {len(failed)}")
    
    # 保存成功模板
    if rows:
        df = pd.DataFrame(rows)
        print(f"\nDataFrame列名: {list(df.columns)}")
        print(f"DataFrame形状: {df.shape}")
        
        # HDF5（aizynthfinder格式）
        hdf5_path = OUTPUT_DIR / 'flavonoid_templates.hdf5'
        df_aizynth = df[['retro_template', 'template_hash', 'classification', 'library_occurence']].copy()
        df_aizynth.to_hdf(hdf5_path, key='table', mode='w', format='table', data_columns=True)
        print(f"\n已保存HDF5: {hdf5_path}")
        
        # 完整CSV
        csv_path = OUTPUT_DIR / 'flavonoid_templates.csv.gz'
        df.to_csv(csv_path, index=False, compression='gzip')
        print(f"已保存CSV: {csv_path}")
        
        # 元数据
        meta_path = OUTPUT_DIR / 'flavonoid_templates_metadata.csv'
        df[['template_code', 'reaction_id', 'name', 'classification', 'priority', 'source_file', 'doi']].to_csv(meta_path, index=False)
        print(f"已保存元数据: {meta_path}")
        
        # 分类统计
        print(f"\n分类分布:\n{df['classification'].value_counts().to_string()}")
    else:
        print("\n警告: 没有成功转换的模板")
    
    # 保存失败报告
    if failed:
        df_failed = pd.DataFrame(failed)
        failed_path = OUTPUT_DIR / 'flavonoid_templates_failed.csv'
        df_failed.to_csv(failed_path, index=False)
        print(f"\n已保存失败报告: {failed_path}")
        print(f"\n前5个失败:")
        for f in failed[:5]:
            print(f"  {f['reaction_id']}: {f['name']} -> {f['error']}")
    
    # 综合报告
    report_path = OUTPUT_DIR / 'conversion_report.txt'
    with open(report_path, 'w') as f:
        f.write(f"总模板: {len(all_templates)}\n")
        f.write(f"成功: {len(rows)}\n")
        f.write(f"失败: {len(failed)}\n")
        if rows:
            f.write(f"\n分类分布:\n{df['classification'].value_counts().to_string()}\n")
    print(f"\n已保存报告: {report_path}")


if __name__ == '__main__':
    main()
