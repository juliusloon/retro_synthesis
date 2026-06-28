#!/usr/bin/env python3
"""
Generate AiZynthFinder-compatible biosynthesis template JSON and stock
from the extracted flavonoid biosynthesis reactions.
"""
import json
import csv
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

ROOT = Path("/home/ljx/retro_synthesis")
OUT_DIR = ROOT / "templates" / "flavonoid_biosynthesis"
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


def inchikey(smi: str) -> str:
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return ""
        return Chem.MolToInchiKey(mol)
    except Exception:
        return ""


def classify_reaction(header: str, raw: str) -> str:
    h = header.lower()
    r = raw
    if "二氢查耳酮" in h or ("[OH-]" in r and "C(=O)C=C" in r):
        return "dihydrochalcone_biosynthesis"
    if "异黄酮" in h or "O=[Fe]" in r or "[CH]" in r or "[Fe]=*" in r or "[O+]=C" in r:
        return "isoflavone_biosynthesis"
    if "其他黄酮" in h or "COC1=CC(=C[C](CO))" in r or "[OH+]" in r or "[H-]" in r:
        return "radical_coupling"
    if "[R]C(=O)CC(=O)" in r or (r.count("CC(=O)") >= 2 and "[R]" in r):
        return "polyketide_biosynthesis"
    if "C=CC(=O)c2" in r and "Oc" in r:
        return "chalcone_cyclization"
    if "O=C1CC" in r and ("O=c1cc" in r or "O=c1c(O)c" in r):
        return "core_oxidation"
    return "biosynthesis_general"


def name_for_reaction(header: str, idx: int, classification: str, raw: str) -> str:
    short = header.split(" / ")[-1] if header else classification.replace("_", " ")
    return f"{short} step {idx}"


def split_multi_product_reactions(reactions):
    """For reactions with multiple products, yield one template per product
    (each product appears as the first product so it can be matched in retro)."""
    out = []
    for r in reactions:
        products = [p.strip() for p in r["products_smiles"].split(".") if p.strip()]
        reactants = r["reactants_smiles"]
        if len(products) <= 1:
            out.append((r, reactants, r["products_smiles"]))
            continue
        for i, prod in enumerate(products):
            others = [p for j, p in enumerate(products) if j != i]
            product_str = ".".join([prod] + others)
            out.append((r, reactants, product_str))
    return out


def clean_smarts(raw: str) -> str:
    """Convert raw reaction SMILES/SMARTS into a forward SMARTS usable by RDKit."""
    smarts = raw.strip()
    # Replace generic R groups with SMARTS wildcard (any atom)
    smarts = smarts.replace("[R]", "[*]").replace("[r]", "[*]")
    # Enzyme/cofactor placeholders: keep as wildcard atoms
    smarts = smarts.replace("O=[Fe]-*", "[*]").replace("O[Fe]*", "[*]").replace("[Fe]=*", "[*]")
    smarts = smarts.replace("CC(=C)=CC*", "CC(=C)CC[*]")
    return smarts


def validate_forward_smarts(smarts: str) -> tuple:
    if ">>" not in smarts:
        return False, "no arrow"
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        if rxn is None:
            return False, "None"
        return True, f"reactants={rxn.GetNumReactantTemplates()}, products={rxn.GetNumProductTemplates()}"
    except Exception as e:
        return False, str(e)[:120]


def generate_templates():
    data = json.loads(REACTIONS_JSON.read_text(encoding="utf-8"))
    templates = []
    skipped = []
    for r, reactants, products in split_multi_product_reactions(data["reactions"]):
        raw_forward = f"{reactants}>>{products}"
        smarts = clean_smarts(raw_forward)
        ok, msg = validate_forward_smarts(smarts)
        if not ok:
            skipped.append({"reaction_id": r["reaction_id"], "error": msg, "smarts": smarts})
            continue
        classification = classify_reaction(r["header"], r["raw_smiles"])
        templates.append({
            "reaction_id": r["reaction_id"],
            "name": name_for_reaction(r["header"], len(templates)+1, classification, r["raw_smiles"]),
            "category": classification,
            "smarts": smarts,
            "conditions": "enzymatic / biosynthetic",
            "stereoselectivity": "retained where specified",
            "substrate_scope": r["header"] or "flavonoid biosynthesis",
            "citation": "From 黄酮化合物的生物合成.md",
            "doi": "",
            "priority": "medium",
            "notes": f"Original: {r['raw_smiles']}; header: {r['header']}",
        })

    # Write single JSON file
    out_json = OUT_DIR / "templates_biosynthesis.json"
    out_json.write_text(json.dumps(templates, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote {len(templates)} biosynthesis templates to {out_json}")
    if skipped:
        print(f"Skipped {len(skipped)} invalid reactions")
        for s in skipped[:5]:
            print(f"  {s['reaction_id']}: {s['error']}")
    return templates


def generate_stock():
    data = json.loads(REACTIONS_JSON.read_text(encoding="utf-8"))
    seen = set()
    rows = []

    def add(smi, name, category, subcategory, notes, price_tier="$"):
        smi = canonical(smi)
        if not smi or smi in seen or "*" in smi or "[R]" in smi:
            return
        seen.add(smi)
        mw = 0.0
        try:
            mw = round(Descriptors.MolWt(Chem.MolFromSmiles(smi)), 2)
        except Exception:
            pass
        rows.append({
            "name": name,
            "smiles": smi,
            "mw": mw,
            "category": category,
            "subcategory": subcategory,
            "notes": notes,
            "price_tier": price_tier,
        })

    # Compounds from the document
    for c in data.get("compounds", []):
        can = c.get("canonical", "")
        if can:
            add(can, "biosynthesis_example", "flavonoid_biosynthesis_example", "example", f"raw: {c.get('raw','')}")

    # Common biosynthetic precursors / cofactors
    cofactors = [
        ("O=C(O)/C=C/c1ccc(O)cc1", "4-coumaric acid", "cofactor", "starter", "phenylpropanoid starter"),
        ("O=C(O)CC(=O)O", "malonic acid", "cofactor", "extender", "polyketide extender surrogate"),
        ("O=C(O)C(=O)CCC(=O)O", "2-oxoglutaric acid", "cofactor", "cosubstrate", "AKG / 2-酮戊二酸"),
        ("O=O", "oxygen", "cofactor", "oxidant", "O2"),
        ("O=C=O", "carbon dioxide", "cofactor", "byproduct", "CO2"),
        ("O", "water", "cofactor", "solvent", "H2O"),
        ("COc1cc(/C=C/CO)ccc1O", "coniferyl alcohol", "cofactor", "phenylpropanoid", "lignin precursor"),
        ("CC(C)=CC=C", "isoprene unit", "cofactor", "prenyl_donor", "DMAPP surrogate"),
        ("Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1", "resveratrol", "flavonoid_biosynthesis_example", "stilbene", "stilbene example"),
        ("O=C1C[C@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21", "naringenin", "flavonoid_biosynthesis_example", "flavanone", "naringenin"),
        ("O=C1c2c(O)cc(O)cc2O[C@H](c2ccc(O)c(O)c2)[C@H]1O", "taxifolin", "flavonoid_biosynthesis_example", "dihydroflavonol", "taxifolin"),
    ]
    for smi, name, cat, subcat, note in cofactors:
        add(smi, name, cat, subcat, note)

    stock_path = OUT_DIR / "flavonoid_biosynthesis_stock.csv"
    with open(stock_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "smiles", "mw", "category", "subcategory", "notes", "price_tier"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} stock entries to {stock_path}")

    keys = sorted(set(k for k in [inchikey(r["smiles"]) for r in rows] if k))
    key_path = OUT_DIR / "flavonoid_biosynthesis_stock_inchikeys.txt"
    key_path.write_text("\n".join(keys) + "\n", encoding="utf-8")
    print(f"Wrote {len(keys)} InChIKeys to {key_path}")


def main():
    generate_templates()
    generate_stock()


if __name__ == "__main__":
    main()
