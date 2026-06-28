#!/usr/bin/env python3
import json
import glob
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
    if not smarts_str or not isinstance(smarts_str, str):
        return False, "empty or not string"
    if '>>' not in smarts_str:
        return False, "no reaction arrow"
    try:
        rxn = Chem.ReactionFromSmarts(smarts_str)
        if rxn is None:
            return False, "ReactionFromSmarts returned None"
        if rxn.GetNumReactantTemplates() == 0:
            return False, "no reactant templates"
        if rxn.GetNumProductTemplates() == 0:
            return False, "no product templates"
        return True, "OK"
    except Exception as e:
        return False, str(e)

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

# 分析每个模板的解析情况
for t in all_templates[:5]:  # 先看前5个
    rid = t.get('reaction_id')
    print(f"=== {rid} ===")
    print(f"name: {t.get('name')}")
    
    smarts = t.get('smarts', '')
    react = t.get('reaction_smarts', '')
    
    ok1, msg1 = parse_smarts(smarts)
    ok2, msg2 = parse_smarts(react)
    
    print(f"smarts 解析: {ok1} ({msg1})")
    print(f"reaction_smarts 解析: {ok2} ({msg2})")
    print(f"smarts 前80: {smarts[:80]}")
    print(f"reaction_smarts 前80: {react[:80]}\n")

# 统计失败原因
failure_reasons = {}
for t in all_templates:
    smarts = t.get('smarts', '')
    react = t.get('reaction_smarts', '')
    ok1, msg1 = parse_smarts(smarts)
    ok2, msg2 = parse_smarts(react)
    if not ok1 and not ok2:
        key = f"smarts:{msg1}|reaction_smarts:{msg2}"
        failure_reasons[key] = failure_reasons.get(key, 0) + 1

print("\n=== 失败原因统计 ===")
for k, v in sorted(failure_reasons.items(), key=lambda x: -x[1])[:10]:
    print(f"{v}: {k}")

print(f"\n总失败数: {sum(failure_reasons.values())}")
