#!/usr/bin/env python3
import json
from rdkit.Chem import AllChem
from pathlib import Path

TEMPLATE_FILES = [
    "templates_glycosylation.json",
    "templates_deprotection.json",
    "templates_cyclization.json",
    "templates_protection_synthesis.json",
    "templates_additional.json",
]

def parse_smarts(smarts_str):
    if not smarts_str or not isinstance(smarts_str, str):
        return False, "empty"
    if '>>' not in smarts_str:
        return False, "no arrow"
    try:
        rxn = AllChem.ReactionFromSmarts(smarts_str)
        if rxn is None:
            return False, "None"
        if rxn.GetNumReactantTemplates() == 0 or rxn.GetNumProductTemplates() == 0:
            return False, "missing templates"
        return True, "OK"
    except Exception as e:
        return False, str(e)[:80]

templates_dir = Path("/home/ljx/retro_synthesis/templates")
all_templates = []

for fname in TEMPLATE_FILES:
    fpath = templates_dir / fname
    with open(fpath, 'r', encoding='utf-8') as f:
        data = json.load(f)
    items = data if isinstance(data, list) else data.get("templates", [])
    for item in items:
        item['source_file'] = fname
    all_templates.extend(items)

print(f"总计: {len(all_templates)} 个模板\n")

success = 0
failed = 0
failed_details = []

for t in all_templates:
    rid = t.get('reaction_id')
    name = t.get('name', '')
    smarts = t.get('smarts', '')
    react = t.get('reaction_smarts', '')
    
    ok1, msg1 = parse_smarts(smarts)
    ok2, msg2 = parse_smarts(react)
    
    if ok1 or ok2:
        success += 1
    else:
        failed += 1
        failed_details.append((rid, name, msg1, msg2))

print(f"成功: {success}")
print(f"失败: {failed}")

if failed > 0:
    print(f"\n前10个失败示例:")
    for rid, name, msg1, msg2 in failed_details[:10]:
        print(f"\n{rid}: {name}")
        print(f"  smarts: {msg1}")
        print(f"  reaction_smarts: {msg2}")
