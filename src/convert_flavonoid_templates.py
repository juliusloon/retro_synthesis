#!/usr/bin/env python3
"""
将同事建立的JSON模板转换为aizynthfinder兼容的HDF5+CSV格式
"""

import json
import glob
import hashlib
import pandas as pd
from rdkit import Chem
from pathlib import Path

TEMPLATE_FILES = [
    "templates_glycosylation.json",
    "templates_deprotection.json",
    "templates_cyclization.json",
    "templates_protection_synthesis.json",
    "templates_additional.json",
]

def parse_smarts(smarts_str):
    """尝试将SMARTS解析为反应，返回是否成功"""
    if not smarts_str or not isinstance(smarts_str, str):
        return False, None
    try:
        rxn = Chem.ReactionFromSmarts(smarts_str)
        if rxn is None:
            return False, None
        # 验证反应有产物和反应物
        if rxn.GetNumReactantTemplates() == 0 or rxn.GetNumProductTemplates() == 0:
            return False, None
        return True, rxn
    except Exception as e:
        return False, None


def load_all_templates():
    """加载并合并所有JSON模板"""
    templates_dir = Path("/home/ljx/retro_synthesis/templates")
    all_templates = []
    
    for fname in TEMPLATE_FILES:
        fpath = templates_dir / fname
        if not fpath.exists():
            print(f"警告: 文件不存在 {fpath}")
            continue
        with open(fpath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # 有些JSON可能是字典（带metadata），有些是列表
        if isinstance(data, dict):
            items = data.get("templates", [])
        elif isinstance(data, list):
            items = data
        else:
            print(f"警告: {fname} 格式未知")
            continue
        
        print(f"{fname}: {len(items)} 个模板")
        for item in items:
            item['source_file'] = fname
        all_templates.extend(items)
    
    print(f"\n总计加载: {len(all_templates)} 个模板")
    return all_templates


def convert_to_aizynth_format(templates):
    """转换为aizynthfinder格式"""
    rows = []
    stats = {'success': 0, 'fallback': 0, 'failed': 0}
    
    for idx, t in enumerate(templates):
        reaction_id = t.get('reaction_id', f'TMP_{idx:03d}')
        name = t.get('name', '')
        category = t.get('category', '')
        
        # 优先使用 smarts，失败则使用 reaction_smarts
        smarts = t.get('smarts', '')
        react_smarts = t.get('reaction_smarts', '')
        
        chosen_field = 'smarts'
        parsed, rxn = parse_smarts(smarts)
        
        if not parsed and react_smarts:
            chosen_field = 'reaction_smarts'
            parsed, rxn = parse_smarts(react_smarts)
            if parsed:
                stats['fallback'] += 1
        else:
            stats['success'] += 1
        
        if not parsed:
            stats['failed'] += 1
            print(f"  解析失败 [{reaction_id}]: {name}")
            print(f"    smarts: {smarts[:80]}...")
            if react_smarts:
                print(f"    reaction_smarts: {react_smarts[:80]}...")
            continue
        
        # 规范化SMARTS字符串
        # RDKit可能重排格式，我们用原始字符串保持原子映射
        final_smarts = t.get(chosen_field, '')
        
        # 生成template_hash
        template_hash = hashlib.sha256(final_smarts.encode()).hexdigest()
        
        # 分类：把category映射到类似USPTO的classification
        classification = category if category else 'flavonoid_custom'
        
        row = {
            'template_code': idx,
            'retro_template': final_smarts,
            'template_hash': template_hash,
            'classification': classification,
            'library_occurence': 1,
            'reaction_id': reaction_id,
            'name': name,
            'source_file': t.get('source_file', ''),
            'priority': t.get('priority', 'medium'),
            'doi': t.get('doi', ''),
            'conditions': t.get('conditions', ''),
            'stereoselectivity': t.get('stereoselectivity', ''),
            'notes': t.get('notes', ''),
        }
        rows.append(row)
    
    df = pd.DataFrame(rows)
    print(f"\n=== 解析统计 ===")
    print(f"成功 (smarts): {stats['success']}")
    print(f"回退 (reaction_smarts): {stats['fallback']}")
    print(f"失败: {stats['failed']}")
    
    return df


def save_outputs(df):
    """保存HDF5和CSV"""
    output_dir = Path("/home/ljx/retro_synthesis/templates")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 保存主HDF5（aizynthfinder使用）
    hdf5_path = output_dir / 'flavonoid_templates.hdf5'
    df_aizynth = df[['retro_template', 'template_hash', 'classification', 'library_occurence']].copy()
    df_aizynth.to_hdf(hdf5_path, key='table', mode='w', format='table', data_columns=True)
    print(f"\n已保存HDF5: {hdf5_path}")
    
    # 保存CSV
    csv_path = output_dir / 'flavonoid_templates.csv.gz'
    df.to_csv(csv_path, index=False, compression='gzip')
    print(f"已保存CSV: {csv_path}")
    
    # 保存元数据CSV（不含大段文本，便于查看）
    meta_path = output_dir / 'flavonoid_templates_metadata.csv'
    meta_cols = ['template_code', 'reaction_id', 'name', 'classification', 'priority', 'source_file', 'doi']
    df[meta_cols].to_csv(meta_path, index=False)
    print(f"已保存元数据: {meta_path}")
    
    # 统计报告
    report_path = output_dir / 'conversion_report.txt'
    with open(report_path, 'w') as f:
        f.write("黄酮类模板转换报告\n")
        f.write("="*60 + "\n")
        f.write(f"总模板数: {len(df)}\n")
        f.write(f"分类分布:\\n")
        f.write(df['classification'].value_counts().to_string())
        f.write(f"\\n\\n优先级分布:\\n")
        f.write(df['priority'].value_counts().to_string())
    print(f"已保存报告: {report_path}")


if __name__ == '__main__':
    templates = load_all_templates()
    df = convert_to_aizynth_format(templates)
    save_outputs(df)
    print(f"\n完成！")
