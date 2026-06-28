#!/usr/bin/env python3
"""
Generate structural-diversity template JSON files and stock library
from the extracted flavonoid structure catalog.
"""
import json
import re
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

ROOT = Path("/home/ljx/retro_synthesis")
OUT_DIR = ROOT / "templates" / "flavonoid_structural_diversity"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CLASS_JSON = OUT_DIR / "flavonoid_structural_classes.json"


def canonical(smi: str) -> str:
    """Canonical SMILES; return empty string if invalid or query."""
    if not smi or "*" in smi or "[R]" in smi:
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


def validate_forward_smarts(smarts: str) -> tuple:
    """Return (ok, error)."""
    if ">>" not in smarts:
        return False, "no arrow"
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        if rxn is None:
            return False, "None"
        return True, f"reactants={rxn.GetNumReactantTemplates()}, products={rxn.GetNumProductTemplates()}"
    except Exception as e:
        return False, str(e)[:120]


# ---------------------------------------------------------------------------
# Template definitions (forward direction: reactants >> products)
# ---------------------------------------------------------------------------

def T(rid, name, category, smarts, conditions="", stereo="", scope="", citation="", doi="", priority="medium", notes=""):
    ok, msg = validate_forward_smarts(smarts)
    if not ok:
        print(f"WARN {rid} invalid SMARTS: {msg}\n  {smarts}")
    return {
        "reaction_id": rid,
        "name": name,
        "category": category,
        "smarts": smarts,
        "conditions": conditions,
        "stereoselectivity": stereo,
        "substrate_scope": scope,
        "citation": citation,
        "doi": doi,
        "priority": priority,
        "notes": notes,
    }


# Shared accurate chromone/isoflavone core SMARTS (atom-mapped)
FLAVONE_CORE_PRODUCT = "[#8:1]=[#6:2]1:[#6:3]:[#6:4](-[#6:5]2:[#6:6]:[#6:7]:[#6:8]:[#6:9]:[#6:10]:2):[#8:11]:[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:1:2"
FLAVANONE_CORE_REACTANT = "[#8:1]=[#6:2]1-[#6:3]-[#6:4](-[#6:5]2:[#6:6]:[#6:7]:[#6:8]:[#6:9]:[#6:10]:2)-[#8:11]-[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:2-1"
FLAVAN_CORE_REACTANT = "[#6:2]1-[#6:3]-[#6:4](-[#6:5]2:[#6:6]:[#6:7]:[#6:8]:[#6:9]:[#6:10]:2)-[#8:11]-[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:2-1"
FLAVONOL_PRODUCT = "[#8:1]=[#6:2]1:[#6:3](-[#8:18]):[#6:4](-[#6:5]2:[#6:6]:[#6:7]:[#6:8]:[#6:9]:[#6:10]:2):[#8:11]:[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:1:2"
DIHYDROFLAVONOL_PRODUCT = "[#8:1]=[#6:2]1-[#6:3](-[#8:18])-[#6:4](-[#6:5]2:[#6:6]:[#6:7]:[#6:8]:[#6:9]:[#6:10]:2)-[#8:11]-[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:2-1"
ISOFLAVONE_REACTANT = "[#8:1]=[#6:2]1:[#6:3](-[#6:4]2:[#6:5]:[#6:6]:[#6:7]:[#6:8]:[#6:9]:2):[#6:10]:[#8:11]:[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:1:2"
ISOFLAVANONE_PRODUCT = "[#8:1]=[#6:2]1-[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:2-[#8:11]-[#6:10]-[#6:3]-1-[#6:4]2:[#6:5]:[#6:6]:[#6:7]:[#6:8]:[#6:9]:2"
ISOFLAVAN_REACTANT = "[#6:2]1-[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:2-[#8:11]-[#6:10]-[#6:3]-1-[#6:4]2:[#6:5]:[#6:6]:[#6:7]:[#6:8]:[#6:9]:2"

STRUCTURAL_CLASS_TEMPLATES = [
    T(
        "SC_001", "Chalcone reduction to dihydrochalcone", "structural_class",
        "[c:1][C:2](=[O:3])[C:4]=[C:5][c:6]>>[c:1][C:2](=[O:3])[C:4][C:5][c:6]",
        conditions="H2, Pd/C, EtOAc or MeOH",
        scope="alpha,beta-unsaturated ketones",
        citation="General hydrogenation",
        priority="high",
        notes="Retro: oxidation of dihydrochalcone back to chalcone",
    ),
    T(
        "SC_002", "Chalcone cyclization to flavanone", "structural_class",
        "[#8:1]=[#6:2](-[#6:3]1:[#6:4]:[#6:5]:[#6:6]:[#6:7]:[#6:8]:1)/[#6:9]=[#6:10]/[#6:11]1:[#6:12]:[#6:13]:[#6:14]:[#6:15]:[#6:16]:1-[#8:17]>>[#8:1]=[#6:2]1-[#6:9]-[#6:10](-[#6:11]2:[#6:12]:[#6:13]:[#6:14]:[#6:15]:[#6:16]:2)-[#8:17]-[#6:3]2:[#6:4]:[#6:5]:[#6:6]:[#6:7]:[#6:8]:2-1",
        conditions="NaOAc, EtOH, reflux or acid/base cyclization",
        scope="2'-hydroxychalcones",
        citation="Claisen-Schmidt cyclization",
        priority="high",
        notes="6-membered C-ring closure; retro opens flavanone to 2'-hydroxychalcone",
    ),
    T(
        "SC_003", "Flavanone dehydrogenation to flavone", "structural_class",
        f"{FLAVANONE_CORE_REACTANT}>>{FLAVONE_CORE_PRODUCT}",
        conditions="DDQ, benzene, reflux or Pd/C dehydrogenation",
        scope="flavanones -> flavones",
        priority="high",
    ),
    T(
        "SC_004", "Flavone C3 hydroxylation to flavonol", "structural_class",
        f"{FLAVONE_CORE_PRODUCT}.[#8:18]>>{FLAVONOL_PRODUCT}",
        conditions="m-CPBA or osmium oxidation / hydrolysis",
        scope="flavones -> flavonols",
        priority="medium",
    ),
    T(
        "SC_005", "Flavanone C3 hydroxylation to dihydroflavonol", "structural_class",
        f"{FLAVANONE_CORE_REACTANT}.[#8:18]>>{DIHYDROFLAVONOL_PRODUCT}",
        conditions="m-CPBA, CH2Cl2, 0°C",
        scope="flavanones -> dihydroflavonols",
        priority="medium",
    ),
    T(
        "SC_006", "Flavan oxidation to flavanone", "structural_class",
        f"{FLAVAN_CORE_REACTANT}.[#8:1]>>{FLAVANONE_CORE_REACTANT}",
        conditions="NaBH4, MeOH, 0°C then reoxidation sequence",
        scope="flavans -> flavanones",
        priority="low",
        notes="Retro deoxygenates flavanone to flavan",
    ),
    T(
        "SC_007", "2'-Hydroxychalcone oxidative cyclization to aurone", "structural_class",
        "[#8:1]=[#6:2](-[#6:3]1:[#6:4]:[#6:5]:[#6:6]:[#6:7]:[#6:8]:1)/[#6:9]=[#6:10]/[#6:11]1:[#6:12]:[#6:13]:[#6:14]:[#6:15]:[#6:16]:1-[#8:17]>>[#8:1]=[#6:2]1-[#6:9](=[#6:10][#6:3]2:[#6:4]:[#6:5]:[#6:6]:[#6:7]:[#6:8]:2)[#6:11]2:[#6:12]:[#6:13]:[#6:14]:[#6:15]:[#6:16]:2-[#8:17]-1",
        conditions="H2O2, NaOH, MeOH",
        scope="2'-hydroxychalcones -> aurones",
        priority="medium",
    ),
    T(
        "SC_008", "2'-Hydroxychalcone isomerization to isoflavone", "structural_class",
        "[#8:1]=[#6:2](-[#6:12]1:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:1)/[#6:10]=[#6:3]/[#6:4]1:[#6:5]:[#6:6]:[#6:7]:[#6:8]:[#6:9]:1-[#8:11]>>[#8:1]=[#6:2]1:[#6:3](-[#6:4]2:[#6:5]:[#6:6]:[#6:7]:[#6:8]:[#6:9]:2):[#6:10]:[#8:11]:[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:1:2",
        conditions="Base-catalyzed rearrangement; Algar-Flynn-Oyamada",
        scope="2'-hydroxychalcones -> isoflavones",
        priority="medium",
    ),
    T(
        "SC_009", "Isoflavone reduction to isoflavanone", "structural_class",
        f"{ISOFLAVONE_REACTANT}>>{ISOFLAVANONE_PRODUCT}",
        conditions="H2, Pd/C or NaBH4",
        scope="isoflavones -> isoflavanones",
        priority="medium",
    ),
    T(
        "SC_010", "Isoflavan oxidation to isoflavanone", "structural_class",
        f"{ISOFLAVAN_REACTANT}.[#8:1]>>{ISOFLAVANONE_PRODUCT}",
        conditions="NaBH4, EtOH then reoxidation",
        scope="isoflavans -> isoflavanones",
        priority="low",
        notes="Retro deoxygenates isoflavanone to isoflavan",
    ),
    T(
        "SC_015", "Isoflavanone oxidation to isoflavone", "structural_class",
        f"{ISOFLAVANONE_PRODUCT}>>{ISOFLAVONE_REACTANT}",
        conditions="DDQ or Pd/C dehydrogenation",
        scope="isoflavanones -> isoflavones",
        priority="medium",
    ),
    T(
        "SC_011", "Isoflavone cyclization to pterocarpan", "structural_class",
        "[#8:1]=[#6:2]1:[#6:3](-[#6:4]2:[#6:5]:[#6:6]:[#6:7]:[#6:8]:[#6:9]:2):[#6:10]:[#8:11]:[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:1:2>>[#6:2]1-[#6:3]2-[#6:4]1:[#6:5]:[#6:6]:[#6:7]:[#6:8]:[#6:9]:2-[#8:11]-[#6:10]-[#6:12]1:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:1",
        conditions="NaBH4 reduction then acid cyclization",
        scope="isoflavones with B-ring o-hydroxy -> pterocarpans",
        priority="low",
    ),
    T(
        "SC_012", "Flavone oxidative cyclization to coumestan", "structural_class",
        "[#8:1]=[#6:2]1:[#6:3]:[#6:4](-[#6:5]2:[#6:6]:[#6:7]:[#6:8]:[#6:9]:[#6:10]:2):[#8:11]:[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:1:2>>[#8:1]=[#6:2]1:[#6:3]:[#6:4](-[#6:5]2:[#6:6]:[#6:7]:[#6:8]:[#6:9]:[#6:10]:2):[#8:11]:[#6:12]2:[#6:13]:[#6:14]:[#6:15]:[#6:16]:[#6:17]:1:2",
        conditions="Oxidative phenolic coupling, DDQ or hypervalent iodine",
        scope="7-hydroxyisoflavones / coumestans",
        priority="low",
        notes="Placeholder: coumestan formation is highly specific; kept for manual curation",
    ),
]

GLYCOSIDE_TEMPLATES = [
    T(
        "GLY_101", "Phenolic O-glycosylation with glycosyl donor", "glycosylation",
        "[c:1][O:2][H:3].[C:4]1[C:5][C:6][C:7][C:8][O:9]1>>[c:1][O:2][C:4]1[C:5][C:6][C:7][C:8][O:9]1",
        conditions="TMSOTf or BF3.OEt2, 4A MS, CH2Cl2, -78°C to 0°C",
        stereo="1,2-trans with participating C2 protecting group",
        scope="flavonoid phenols + protected glycosyl donor",
        citation="Schmidt glycosylation",
        doi="10.1002/anie.198602121",
        priority="high",
        notes="Generic donor represented as a pyranose ring; leaving group omitted for matching",
    ),
    T(
        "GLY_102", "Disaccharide O-glycosylation (rutinoside type)", "glycosylation",
        "[c:1][O:2][H:3].[C:4]1[C:5][C:6][C:7][C:8]([O:9][C:10]2[C:11][C:12][C:13][C:14][O:9]2)[O:15]1>>[c:1][O:2][C:4]1[C:5][C:6][C:7][C:8]([O:9][C:10]2[C:11][C:12][C:13][C:14][O:9]2)[O:15]1",
        conditions="TMSOTf, CH2Cl2, -40°C",
        scope="flavonoid 7-OH + disaccharide donor",
        priority="high",
        notes="1,6-linked disaccharide donor generic representation",
    ),
    T(
        "GLY_103", "C-glycosylation of flavonoid A-ring", "c_glycosylation",
        "[c:1]1[c:2][c:3]([O:4])[c:5][c:6]1.[C:7]1[C:8][C:9][C:10][C:11][O:12]1>>[c:1]1[c:2][c:3]([O:4])[c:5]([C:7]2[C:8][C:9][C:10][C:11][O:12]2)[c:6]1",
        conditions="Lewis acid, CH2Cl2, -78°C to rt",
        scope="electron-rich flavonoid A-ring C6/C8 + glycosyl donor",
        priority="medium",
        notes="C-C glycosidic bond formation",
    ),
]

SUBSTITUENT_TEMPLATES = [
    T(
        "SUB_101", "Phenol O-methylation", "substitution",
        "[c:1][O:2][H:3].[C:4][I:5]>>[c:1][O:2][C:4].[I:5]",
        conditions="K2CO3, MeI, acetone, reflux",
        scope="phenols",
        priority="high",
    ),
    T(
        "SUB_102", "Phenol O-acetylation", "substitution",
        "[c:1][O:2][H:3].[C:4](=[O:5])[Cl:6]>>[c:1][O:2][C:4](=[O:5]).[Cl:6]",
        conditions="Et3N, AcCl or Ac2O, CH2Cl2, 0°C to rt",
        scope="phenols",
        priority="high",
    ),
    T(
        "SUB_103", "Phenol O-benzoylation", "substitution",
        "[c:1][O:2][H:3].[C:4](=[O:5])([Cl:6])[c:7]1ccccc1>>[c:1][O:2][C:4](=[O:5])[c:7]1ccccc1.[Cl:6]",
        conditions="Et3N, BzCl, CH2Cl2, 0°C to rt",
        scope="phenols",
        priority="medium",
    ),
    T(
        "SUB_104", "Phenol O-sulfation", "substitution",
        "[c:1][O:2][H:3].[S:4](=[O:5])(=[O:6])([O:7])>>[c:1][O:2][S:4](=[O:5])(=[O:6])[O:7]",
        conditions="SO3.pyridine, DMF, 0°C to rt",
        scope="phenols -> aryl sulfates",
        priority="medium",
    ),
    T(
        "SUB_105", "C-prenylation of arene", "substitution",
        "[c:1][H:2].[C:3](=[C:4][C:5])([C:6])[Cl:7]>>[c:1][C:3](=[C:4][C:5])([C:6]).[Cl:7]",
        conditions="Lewis acid, CH2Cl2, rt",
        scope="electron-rich arenes",
        priority="medium",
        notes="Generic prenyl chloride donor (3-methyl-2-butenyl)",
    ),
    T(
        "SUB_106", "O-prenylation of phenol", "substitution",
        "[c:1][O:2][H:3].[C:4](=[C:5][C:6])([C:7])[Cl:8]>>[c:1][O:2][C:4](=[C:5][C:6])([C:7]).[Cl:8]",
        conditions="K2CO3, prenyl bromide/chloride, acetone, reflux",
        scope="phenols",
        priority="medium",
    ),
    T(
        "SUB_107", "C-methylation of arene", "substitution",
        "[c:1][H:2].[C:3][I:4]>>[c:1][C:3].[I:4]",
        conditions="Base or Friedel-Crafts methylation",
        scope="activated arenes",
        priority="low",
    ),
    T(
        "SUB_108", "C-formylation of arene", "substitution",
        "[c:1][H:2].[C:3](=[O:4])[Cl:5]>>[c:1][C:3](=[O:4]).[Cl:5]",
        conditions="Vilsmeier-Haack or Reimer-Tiemann",
        scope="activated arenes",
        priority="low",
    ),
    T(
        "SUB_109", "C-benzyl installation on arene", "substitution",
        "[c:1][H:2].[C:3]([Cl:4])[c:5]1ccccc1>>[c:1][C:3][c:5]1ccccc1.[Cl:4]",
        conditions="Friedel-Crafts benzylation",
        scope="activated arenes",
        priority="low",
    ),
    T(
        "SUB_110", "Methylenedioxy formation from catechol", "substitution",
        "[c:1]1[c:2]([O:3])[c:4]([O:5])[c:6][c:7][c:8]1.[C:9]([Br:10])[Br:11]>>[c:1]1[c:2]([O:3][C:9][O:5]2)[c:4]2[c:6][c:7][c:8]1.[Br:10][Br:11]",
        conditions="CH2Br2, K2CO3, acetone, reflux",
        scope="ortho-dihydroxyarenes",
        priority="low",
    ),
]

PRENYL_CYCLIZATION_TEMPLATES = [
    T(
        "PRE_101", "Prenyl ether 5-membered ring cyclization", "prenyl_cyclization",
        "[c:1][O:2][C:3](=[C:4][C:5])([C:6])[C:7][C:8]([C:9])=C>>[c:1][O:2][C:3]1[C:4]([C:5])([C:6])[C:7][C:8]1[C:9]",
        conditions="Acid-promoted cyclization",
        scope="O-prenyl flavonoids",
        priority="medium",
    ),
    T(
        "PRE_102", "Prenyl ether 6-membered ring cyclization", "prenyl_cyclization",
        "[c:1][O:2][C:3](=[C:4][C:5])([C:6])[C:7][C:8][C:9]([C:10])=C>>[c:1][O:2][C:3]1[C:4]([C:5])([C:6])[C:7][C:8][C:9]1[C:10]",
        conditions="Acid-promoted cyclization",
        scope="O-prenyl flavonoids with extended side chain",
        priority="medium",
    ),
    T(
        "PRE_103", "C-prenyl cyclization to dihydrobenzofuran", "prenyl_cyclization",
        "[c:1]([C:2](=[C:3][C:4])([C:5])[C:6][C:7]([C:8])=C)[c:9]1[c:10]([O:11][H:12])[c:13]cc[c:9]1>>[c:1]1[c:9]2[c:10]([O:11][C:2]3[C:3]([C:4])([C:5])[C:6][C:7]3[C:8])[c:13]cc[c:9]2[c:1]1",
        conditions="Acid or oxidative cyclization",
        scope="C-prenyl flavonoids with adjacent phenol",
        priority="low",
    ),
]

QUINONE_DA_TEMPLATES = [
    T(
        "QD_101", "Catechol oxidation to ortho-quinone", "quinone",
        "[c:1]1[c:2]([O:3][H:4])[c:5]([O:6][H:7])[c:8][c:9]1>>[O:3]=[C:2]1[C:5](=[O:6])[C:8]=[C:9][C:1]=1",
        conditions="NaIO4, CH2Cl2/H2O or Ag2O",
        scope="catechols",
        priority="medium",
    ),
    T(
        "QD_102", "Hydroquinone oxidation to para-quinone", "quinone",
        "[c:1]1[c:2]([O:3][H:4])cc[c:5]([O:6][H:7])[c:2]1>>[O:3]=[C:2]1C=CC(=O)C=C1",
        conditions="Ag2O, benzene or Fremy's salt",
        scope="hydroquinones",
        priority="medium",
    ),
    T(
        "QD_103", "Retro Diels-Alder fragmentation of chalcone-diene adduct", "retro_diels_alder",
        "[C:1]1[C:2]=[C:3][C:4]([C:5](=[O:6])[c:7]1)[C:8]=[C:9][c:10]1ccccc1>>[C:1]=[C:2][C:3]=[C:4].[O:5]=[C:6][c:7]1ccccc1",
        conditions="Thermal retro-Diels-Alder",
        scope="cyclohexene adducts of chalcones",
        priority="medium",
    ),
]

ALL_TEMPLATE_FILES = {
    "templates_structural_class.json": STRUCTURAL_CLASS_TEMPLATES,
    "templates_glycosides.json": GLYCOSIDE_TEMPLATES,
    "templates_substituents.json": SUBSTITUENT_TEMPLATES,
    "templates_prenyl_cyclizations.json": PRENYL_CYCLIZATION_TEMPLATES,
    "templates_quinone_da.json": QUINONE_DA_TEMPLATES,
}


# ---------------------------------------------------------------------------
# Stock generation
# ---------------------------------------------------------------------------

def build_stock():
    catalog = json.loads(CLASS_JSON.read_text(encoding="utf-8"))
    rows = []
    seen = set()

    def add(smi, name, category, subcategory, notes, price_tier="$"):
        smi = canonical(smi)
        if not smi or smi in seen:
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

    # Representative molecules from the document
    for heading, items in catalog["hierarchy"].items():
        short = heading.split(" / ")[-1]
        for item in items:
            smi = item["canonical_smiles"]
            if not smi or "*" in smi or "[R]" in smi:
                continue
            # Avoid generic sugars represented with * or R groups
            add(smi, f"Structure_example_line_{item['line']}", "flavonoid_structural_example", short, f"from {heading}")

    # Additional building blocks / donors / reagents
    building_blocks = [
        ("CC(C)=CCBr", "prenyl bromide", "prenyl_donor", "prenyl_donor", "C5 prenyl donor"),
        ("CC(C)=CCCl", "prenyl chloride", "prenyl_donor", "prenyl_donor", "C5 prenyl donor"),
        ("CC(C)=CCO", "prenyl alcohol", "prenyl_donor", "prenyl_donor", "C5 prenyl alcohol"),
        ("CC(C)=CCCC(C)=CCBr", "geranyl bromide", "prenyl_donor", "prenyl_donor", "C10 geranyl donor"),
        ("CC(=O)Cl", "acetyl chloride", "acyl_donor", "acyl_donor", "acetylation reagent"),
        ("CC(=O)OC(C)=O", "acetic anhydride", "acyl_donor", "acyl_donor", "acetylation reagent"),
        ("O=C(Cl)c1ccccc1", "benzoyl chloride", "acyl_donor", "acyl_donor", "benzoylation reagent"),
        ("O=C(Cl)/C=C/c1ccccc1", "cinnamoyl chloride", "acyl_donor", "acyl_donor", "cinnamoylation reagent"),
        ("O=C(Cl)/C=C/c1ccc(O)cc1", "p-coumaroyl chloride", "acyl_donor", "acyl_donor", "p-coumaroylation reagent"),
        ("O=S(=O)(O)O", "sulfuric acid", "sulfation_reagent", "sulfation_reagent", "sulfation reagent"),
        ("COS(=O)(=O)O", "methyl sulfate", "sulfation_reagent", "sulfation_reagent", "sulfation reagent"),
        ("CI", "methyl iodide", "methylation_reagent", "methylation_reagent", "methylation reagent"),
        ("CBr", "methyl bromide", "methylation_reagent", "methylation_reagent", "methylation reagent"),
        ("ClCc1ccccc1", "benzyl chloride", "benzylation_reagent", "benzylation_reagent", "benzyl donor"),
        ("C1COCCO1", "1,2-dibromoethane surrogate", "methylenedioxy_reagent", "methylenedioxy_reagent", "methylenedioxy formation"),
        ("O=Cc1ccccc1O", "salicylaldehyde", "flavonoid_precursor", "A_ring", "aurone / benzofuran precursor"),
        ("CC(=O)c1ccccc1O", "2-hydroxyacetophenone", "flavonoid_precursor", "A_ring", "flavonoid A-ring precursor"),
        ("Oc1cc(O)cc(O)c1", "phloroglucinol", "flavonoid_precursor", "A_ring", "1,3,5-trihydroxybenzene"),
        ("Oc1cccc(O)c1", "resorcinol", "flavonoid_precursor", "A_ring", "1,3-dihydroxybenzene"),
        ("Oc1ccccc1O", "catechol", "flavonoid_precursor", "A_ring", "1,2-dihydroxybenzene"),
        ("O=Cc1ccc(O)cc1", "4-hydroxybenzaldehyde", "flavonoid_precursor", "B_ring", "B-ring precursor"),
        ("O=Cc1ccc(O)c(O)c1", "protocatechuic aldehyde", "flavonoid_precursor", "B_ring", "catechol B-ring"),
        ("O=C(/C=C/c1ccccc1)c1ccccc1", "chalcone", "flavonoid_precursor", "chalcone", "parent chalcone"),
        ("O=C1CC(c2ccccc2)Oc2ccccc21", "flavanone", "flavonoid_precursor", "flavanone", "parent flavanone"),
        ("O=c1cc(-c2ccccc2)oc2ccccc12", "flavone", "flavonoid_precursor", "flavone", "parent flavone"),
        ("O=c1c(O)c(-c2ccccc2)oc2ccccc12", "flavonol", "flavonoid_precursor", "flavonol", "parent flavonol"),
        ("O=c1c(-c2ccccc2)coc2ccccc12", "isoflavone", "flavonoid_precursor", "isoflavone", "parent isoflavone"),
        ("O=C1c2ccccc2OCC1c1ccccc1", "isoflavanone", "flavonoid_precursor", "isoflavanone", "parent isoflavanone"),
        ("O=C1C(=Cc2ccccc2)Oc2ccccc21", "aurone", "flavonoid_precursor", "aurone", "parent aurone"),
        ("O=c1oc2ccccc2cc1-c1ccccc1O", "3-arylcoumarin", "flavonoid_precursor", "arylcoumarin", "parent 3-arylcoumarin"),
        ("O=c1oc2ccccc2c2oc3ccccc3c12", "coumestan", "flavonoid_precursor", "coumestan", "parent coumestan"),
        ("O=C1c2ccccc2OC(c2ccccc2)C1O", "dihydroflavonol", "flavonoid_precursor", "dihydroflavonol", "parent dihydroflavonol"),
        ("Oc1cc2c(=O)c(-c3ccc(O)cc3)coc2cc1O", "genistein", "flavonoid_precursor", "isoflavone", "4',5,7-trihydroxyisoflavone"),
        ("COc1ccc(C2CC(=O)c3c(O)cc(O)cc3O2)cc1O", "hesperetin", "flavonoid_precursor", "flavanone", "hesperidin aglycone"),
        ("O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)ccc12", "quercetin", "flavonoid_precursor", "flavonol", "quercetin"),
        ("O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)ccc12", "kaempferol", "flavonoid_precursor", "flavonol", "kaempferol"),
    ]

    for smi, name, cat, subcat, note in building_blocks:
        add(smi, name, cat, subcat, note)

    # Sugars for O/C glycosides
    sugars = [
        ("OC1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O", "D-glucose"),
        ("OC1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O", "gentiobiose"),
        ("C[C@@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2CO)[C@H](O)[C@H](O)[C@H]1O", "rutinose"),
        ("C[C@@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O", "L-rhamnose"),
    ]
    for smi, name in sugars:
        add(smi, name, "sugar", "sugar", "glycoside building block")

    return rows


def main():
    # Write template JSON files
    for fname, templates in ALL_TEMPLATE_FILES.items():
        path = OUT_DIR / fname
        path.write_text(json.dumps(templates, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"Wrote {len(templates)} templates to {path}")

    # Write stock
    stock_rows = build_stock()
    import csv
    stock_path = OUT_DIR / "flavonoid_structural_stock.csv"
    with open(stock_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "smiles", "mw", "category", "subcategory", "notes", "price_tier"])
        writer.writeheader()
        writer.writerows(stock_rows)
    print(f"Wrote {len(stock_rows)} stock entries to {stock_path}")

    # InChIKey file
    keys = [inchikey(r["smiles"]) for r in stock_rows if r["smiles"]]
    keys = sorted(set(k for k in keys if k))
    key_path = OUT_DIR / "flavonoid_structural_stock_inchikeys.txt"
    key_path.write_text("\n".join(keys) + "\n", encoding="utf-8")
    print(f"Wrote {len(keys)} InChIKeys to {key_path}")


if __name__ == "__main__":
    main()
