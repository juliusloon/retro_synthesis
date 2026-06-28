#!/usr/bin/env python3
"""
Extract SMILES examples and section context from
references/黄酮化合物的结构.md for use in the AiZynthFinder flavonoid
structural-diversity library.
"""
import json
import re
from pathlib import Path
from collections import OrderedDict
from rdkit import Chem

MD_PATH = Path("/home/ljx/references/黄酮化合物的结构.md")
OUT_DIR = Path("/home/ljx/retro_synthesis/templates/flavonoid_structural_diversity")
OUT_JSON = OUT_DIR / "flavonoid_structural_classes.json"


def canonical_smiles(smi: str) -> str:
    """Return canonical SMILES if RDKit can parse it; otherwise return original."""
    if not smi:
        return ""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return smi
    try:
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return smi


def extract():
    text = MD_PATH.read_text(encoding="utf-8")
    lines = text.splitlines()

    records = []
    current_headings = []

    for lineno, line in enumerate(lines, start=1):
        # Track headings (#, ##, ###, ####)
        m = re.match(r"^(#{1,4})\s+(.+)$", line)
        if m:
            level = len(m.group(1))
            title = m.group(2).strip()
            # Reset lower/equal-level headings
            current_headings = [(l, t) for l, t in current_headings if l < level]
            current_headings.append((level, title))
            continue

        # Find all $smiles=... tokens on the line
        for match in re.finditer(r"\$smiles=([^`\s]+)", line):
            smi_raw = match.group(1)
            # Strip stray trailing punctuation unlikely in SMILES
            smi_raw = re.sub(r"[.,;:!?]+$", "", smi_raw)
            records.append({
                "line": lineno,
                "headings": [t for _, t in current_headings],
                "raw_smiles": smi_raw,
                "canonical_smiles": canonical_smiles(smi_raw),
            })

    # Build a simple hierarchy summary
    hierarchy = OrderedDict()
    for rec in records:
        key = " / ".join(rec["headings"]) if rec["headings"] else "未分类"
        hierarchy.setdefault(key, []).append(rec)

    output = {
        "source": str(MD_PATH),
        "total_smiles": len(records),
        "hierarchy": {
            k: [
                {
                    "line": r["line"],
                    "raw_smiles": r["raw_smiles"],
                    "canonical_smiles": r["canonical_smiles"],
                }
                for r in v
            ]
            for k, v in hierarchy.items()
        },
    }

    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Extracted {len(records)} SMILES into {OUT_JSON}")
    print(f"Top-level categories: {len(hierarchy)}")
    for k, v in list(hierarchy.items())[:10]:
        print(f"  {k}: {len(v)} example(s)")


if __name__ == "__main__":
    extract()
