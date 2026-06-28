#!/usr/bin/env python3
"""
Extract flavonoid biosynthesis reactions and compound SMILES from
references/黄酮化合物的生物合成.md.
"""
import json
import re
from pathlib import Path
from rdkit import Chem

MD_PATH = Path("/home/ljx/references/黄酮化合物的生物合成.md")
OUT_DIR = Path("/home/ljx/retro_synthesis/templates/flavonoid_biosynthesis")
OUT_DIR.mkdir(parents=True, exist_ok=True)
REACTIONS_JSON = OUT_DIR / "flavonoid_biosynthesis_reactions.json"


def canonical(smi: str) -> str:
    if not smi or smi.lower() in ("nan", "none"):
        return ""
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return ""
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return ""


def split_reaction(smi: str):
    """Split reaction SMILES into left and right parts."""
    if ">>" not in smi:
        return None
    left, right = smi.split(">>", 1)
    return left.strip(), right.strip()


def extract():
    text = MD_PATH.read_text(encoding="utf-8")
    lines = text.splitlines()
    reactions = []
    compounds = []
    seen_reactions = set()
    seen_compounds = set()
    current_headings = []
    in_smiles_block = False
    block_lines = []
    block_header = ""

    def record_reaction(raw, source, header):
        if raw in seen_reactions:
            return
        seen_reactions.add(raw)
        parts = split_reaction(raw)
        if not parts:
            return
        left, right = parts
        reactions.append({
            "reaction_id": f"BIO_{len(reactions)+1:03d}",
            "source": source,
            "header": " / ".join([t for _, t in current_headings]) if current_headings else header,
            "raw_smiles": raw,
            "reactants_smiles": left,
            "products_smiles": right,
            "reactants_canonical": [canonical(s) for s in left.split(".") if s.strip()],
            "products_canonical": [canonical(s) for s in right.split(".") if s.strip()],
        })

    def record_compound(smi):
        smi = smi.strip()
        if not smi or ">>" in smi or smi in seen_compounds:
            return
        seen_compounds.add(smi)
        can = canonical(smi)
        if can:
            compounds.append({"raw": smi, "canonical": can})

    for lineno, line in enumerate(lines, start=1):
        # heading tracking
        m = re.match(r"^(#{1,4})\s+(.+)$", line)
        if m:
            level = len(m.group(1))
            title = m.group(2).strip()
            current_headings = [(l, t) for l, t in current_headings if l < level]
            current_headings.append((level, title))
            continue

        # code block with smiles
        if line.strip().startswith("```smiles"):
            in_smiles_block = True
            block_lines = []
            block_header = " / ".join([t for _, t in current_headings])
            continue
        if in_smiles_block and line.strip().startswith("```"):
            in_smiles_block = False
            block_text = " ".join(block_lines).strip()
            if block_text:
                record_reaction(block_text, f"code_block_line_{lineno}", block_header)
            continue
        if in_smiles_block:
            block_lines.append(line.strip())
            continue

        # inline $smiles=...
        for match in re.finditer(r"\$smiles=([^`\s]+)", line):
            raw = match.group(1)
            raw = re.sub(r"[.,;:!?]+$", "", raw)
            if ">>" in raw:
                record_reaction(raw, f"inline_line_{lineno}", "")
            else:
                record_compound(raw)

    # Also treat compounds extracted from reactions as stock
    for r in reactions:
        for s in r["reactants_canonical"] + r["products_canonical"]:
            if s and s not in seen_compounds:
                seen_compounds.add(s)
                compounds.append({"raw": s, "canonical": s})

    output = {
        "source": str(MD_PATH),
        "total_reactions": len(reactions),
        "total_compounds": len(compounds),
        "reactions": reactions,
        "compounds": compounds,
    }
    REACTIONS_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Extracted {len(reactions)} reactions and {len(compounds)} compounds into {REACTIONS_JSON}")


if __name__ == "__main__":
    extract()
