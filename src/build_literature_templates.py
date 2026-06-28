#!/usr/bin/env python3
"""
Build curated flavonoid reaction templates and stock candidates from markdown notes.

The markdown files mix explicit reaction SMILES with prose-only conditions. This
script keeps the rich provenance in JSON and exports the valid subset in the
small AiZynthFinder HDF5 table expected by custom_expansion.UniformTemplateExpansion.
"""
from __future__ import annotations

import csv
import copy
import hashlib
import itertools
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors


ROOT = Path("/home/ljx/retro_synthesis")
SOURCE_DIR = ROOT / "files"
TEMPLATE_DIR = ROOT / "templates"
OUT_DIR = TEMPLATE_DIR / "literature_curated"
BIO_DIR = TEMPLATE_DIR / "flavonoid_biosynthesis"
STRUCTURAL_DIR = TEMPLATE_DIR / "flavonoid_structural_diversity"

RDLogger.DisableLog("rdApp.warning")
RDLogger.DisableLog("rdApp.error")


CONDITION_KEYWORDS = (
    "反应条件",
    "条件",
    "催化",
    "催化剂",
    "产率",
    "收率",
    "温度",
    "室温",
    "低温",
    "回流",
    "加热",
    "溶剂",
    "脱保护",
    "氧化",
    "还原",
    "选择性",
    "酸",
    "碱",
)

SOLVENT_PATTERNS = {
    "DMSO": r"\bDMSO\b",
    "DMF": r"\bDMF\b",
    "DCM/CH2Cl2": r"\b(?:DCM|CH2Cl2|CHCl3)\b|二氯甲烷|氯仿",
    "MeCN": r"\b(?:MeCN|CH3CN)\b|乙腈",
    "THF": r"\bTHF\b|四氢呋喃",
    "toluene": r"\btoluene\b|甲苯",
    "pyridine": r"\b(?:Py|pyridine)\b|吡啶",
    "MeOH": r"\bMeOH\b|甲醇",
    "EtOH": r"\bEtOH\b|乙醇",
    "AcOH": r"\bAcOH\b|醋酸|乙酸",
    "H2O": r"\bH2O\b|水",
}

CATALYST_PATTERNS = {
    "Ag2O": r"\bAg2O\b|氧化银",
    "Ag2CO3": r"\bAg2CO3\b|碳酸银",
    "AgOTf": r"\bAgOTf\b",
    "AgClO4": r"\bAgClO4\b",
    "BF3.Et2O": r"\bBF3(?:·|\.|-)Et2O\b|\bBF3\b",
    "SnCl4": r"\bSnCl4\b",
    "ZnCl2": r"\bZnCl2\b",
    "TMSOTf": r"\bTMSOTf\b",
    "FeCl3": r"\bFeCl3\b",
    "AlCl3": r"\bAlCl3\b",
    "TBAB": r"\bTBAB\b",
    "DDQ": r"\bDDQ\b",
    "Pd/C": r"\bPd/C\b",
    "K2CO3": r"\bK2CO3\b",
    "NaH": r"\bNaH\b|\[Na\]\[H\]",
    "I2": r"\bI2\b",
    "Ti(NO3)3": r"\bTi\(NO3\)3\b",
    "Cp2HfCl2": r"\bCp2HfCl2\b",
    "Cp2ZrCl2": r"\bCp2ZrCl2\b",
}

SUBSTITUENT_LIBRARY: dict[str, dict[str, Any]] = {
    "H": {"type": "no_substituent", "attachment_smarts": "[H]", "smiles": "[H]"},
    "Me": {"type": "alkyl", "attachment_smarts": "[*:1][CH3]", "smiles": "[*:1]C"},
    "alkyl": {"type": "alkyl_class", "attachment_smarts": "[*:1][C;!$(C=O)]", "smiles": ""},
    "OMe": {"type": "alkoxy", "attachment_smarts": "[*:1][O][CH3]", "smiles": "[*:1]OC"},
    "MeO": {"type": "alkoxy", "attachment_smarts": "[*:1][O][CH3]", "smiles": "[*:1]OC"},
    "OH": {"type": "hydroxy", "attachment_smarts": "[*:1][OH]", "smiles": "[*:1]O"},
    "F": {"type": "halogen", "attachment_smarts": "[*:1][F]", "smiles": "[*:1]F"},
    "Cl": {"type": "halogen", "attachment_smarts": "[*:1][Cl]", "smiles": "[*:1]Cl"},
    "Br": {"type": "halogen", "attachment_smarts": "[*:1][Br]", "smiles": "[*:1]Br"},
    "I": {"type": "halogen", "attachment_smarts": "[*:1][I]", "smiles": "[*:1]I"},
    "halogen": {"type": "halogen_class", "attachment_smarts": "[*:1][F,Cl,Br,I]", "smiles": ""},
    "NO2": {"type": "nitro", "attachment_smarts": "[*:1][N+](=O)[O-]", "smiles": "[*:1][N+](=O)[O-]"},
    "CF3": {"type": "trifluoromethyl", "attachment_smarts": "[*:1][C](F)(F)F", "smiles": "[*:1]C(F)(F)F"},
    "aryl": {"type": "aryl_class", "attachment_smarts": "[*:1][c]", "smiles": ""},
    "aryl_ring_substituent": {"type": "aryl_ring_substituent_class", "attachment_smarts": "[*:1][*]", "smiles": ""},
    "phenyl": {"type": "phenyl", "attachment_smarts": "[*:1]c1ccccc1", "smiles": "[*:1]c1ccccc1"},
    "SCoA": {
        "type": "coenzyme_A_thioester",
        "attachment_smarts": "[*:1][S;X2]",
        "smiles": "",
        "note": "Context defines R = -SCoA; represented structurally as the thioester sulfur attachment, not full CoA.",
    },
    "Enz": {"type": "enzyme_ligand", "attachment_smarts": "[*:1][*]", "smiles": "", "exclude_from_ready": True},
    "OPP": {
        "type": "diphosphate_leaving_group",
        "attachment_smarts": "[*:1][O][P](=O)([O-])[O][P](=O)([O-])[O-]",
        "smiles": "[*:1]OP(=O)([O-])OP(=O)([O-])[O-]",
    },
    "Glu": {
        "type": "beta_D_glucopyranosyl",
        "attachment_smarts": "[*:1][C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O",
        "smiles": "[*:1][C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O",
    },
    "Glc": {
        "type": "beta_D_glucopyranosyl",
        "attachment_smarts": "[*:1][C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O",
        "smiles": "[*:1][C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O",
    },
    "Rha": {
        "type": "alpha_L_rhamnopyranosyl",
        "attachment_smarts": "[*:1][C@H]1O[C@H](C)[C@@H](O)[C@H](O)[C@@H]1O",
        "smiles": "[*:1][C@H]1O[C@H](C)[C@@H](O)[C@H](O)[C@@H]1O",
    },
    "Api": {
        "type": "apiofuranosyl",
        "attachment_smarts": "[*:1][C@H]1OC(CO)(CO)[C@H]1O",
        "smiles": "[*:1][C@H]1OC(CO)(CO)[C@H]1O",
    },
    "Xyl": {
        "type": "xylopyranosyl",
        "attachment_smarts": "[*:1][C@H]1OC[C@H](O)[C@@H](O)[C@H]1O",
        "smiles": "[*:1][C@H]1OC[C@H](O)[C@@H](O)[C@H]1O",
    },
    "Gal": {
        "type": "galactopyranosyl",
        "attachment_smarts": "[*:1][C@H]1O[C@@H](CO)[C@H](O)[C@H](O)[C@H]1O",
        "smiles": "[*:1][C@H]1O[C@@H](CO)[C@H](O)[C@H](O)[C@H]1O",
    },
    "glucuronide": {
        "type": "glucuronopyranosyl",
        "attachment_smarts": "[*:1][C@H]1O[C@@H](C(=O)O)[C@H](O)[C@@H](O)[C@H]1O",
        "smiles": "[*:1][C@H]1O[C@@H](C(=O)O)[C@H](O)[C@@H](O)[C@H]1O",
    },
    "Ac": {"type": "acetyl", "attachment_smarts": "[*:1]C(=O)C", "smiles": "[*:1]C(=O)C"},
    "malonyl": {"type": "malonyl", "attachment_smarts": "[*:1]C(=O)CC(=O)O", "smiles": "[*:1]C(=O)CC(=O)O"},
    "caffeoyl": {
        "type": "caffeoyl",
        "attachment_smarts": "[*:1]C(=O)/C=C/c1ccc(O)c(O)c1",
        "smiles": "[*:1]C(=O)/C=C/c1ccc(O)c(O)c1",
    },
    "cinnamoyl": {
        "type": "cinnamoyl",
        "attachment_smarts": "[*:1]C(=O)/C=C/c1ccccc1",
        "smiles": "[*:1]C(=O)/C=C/c1ccccc1",
    },
    "p-coumaroyl": {
        "type": "p_coumaroyl",
        "attachment_smarts": "[*:1]C(=O)/C=C/c1ccc(O)cc1",
        "smiles": "[*:1]C(=O)/C=C/c1ccc(O)cc1",
    },
    "trimethylsilyl": {"type": "silyl", "attachment_smarts": "[*:1][Si](C)(C)C", "smiles": "[*:1][Si](C)(C)C"},
    "hexyl": {"type": "alkyl", "attachment_smarts": "[*:1]CCCCCC", "smiles": "[*:1]CCCCCC"},
    "octyl": {"type": "alkyl", "attachment_smarts": "[*:1]CCCCCCCC", "smiles": "[*:1]CCCCCCCC"},
    "allyl": {"type": "alkenyl", "attachment_smarts": "[*:1]CC=C", "smiles": "[*:1]CC=C"},
    "butenyl": {"type": "alkenyl", "attachment_smarts": "[*:1]CC=CC", "smiles": "[*:1]CC=CC"},
    "styryl": {"type": "alkenyl_aryl", "attachment_smarts": "[*:1]C=Cc1ccccc1", "smiles": "[*:1]C=Cc1ccccc1"},
}


ATTACHMENT_FRAGMENTS: dict[str, str] = {
    "H": "[H]",
    "Me": "*C",
    "Et": "*CC",
    "i-Pr": "*C(C)C",
    "OH": "*O",
    "OMe": "*OC",
    "MeO": "*OC",
    "F": "*F",
    "Cl": "*Cl",
    "Br": "*Br",
    "I": "*I",
    "NO2": "*[N+](=O)[O-]",
    "CF3": "*C(F)(F)F",
    "SCoA": "*S",
    "OPP": "*OP(=O)([O-])OP(=O)([O-])[O-]",
    "TMS": "*[Si](C)(C)C",
    "hexyl": "*CCCCCC",
    "octyl": "*CCCCCCCC",
    "allyl": "*CC=C",
    "butenyl": "*CC=CC",
    "styryl": "*C=Cc1ccccc1",
    "butadienyl": "*/C=C/C=C",
    "Glc": "*[C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O",
    "Glu": "*[C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O",
    "Rha": "*[C@H]1O[C@H](C)[C@@H](O)[C@H](O)[C@@H]1O",
    "Api": "*[C@H]1OC(CO)(CO)[C@H]1O",
    # Neohesperidosyl-like Glc-Rha. The source text uses Glc-Rha as a shorthand;
    # this encodes an O-linked rhamnosyl substituent on the glucose core.
    "Glc-Rha": "*[C@H]1O[C@@H](CO)[C@H](O)[C@@H](O[C@H]2O[C@H](C)[C@@H](O)[C@H](O)[C@@H]2O)[C@H]1O",
    # Kaempferol glycoside shorthand from the notes: glucose attached to the
    # aglycone, one extra glucose, and a caffeoyl ester on the primary glucose.
    "Glc(3-Glc)-4-caffeoyl": "*[C@H]1O[C@@H](CO)[C@H](O[C@H]2O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]2O)[C@@H](OC(=O)/C=C/c2ccc(O)c(O)c2)[C@H]1O",
    "Glc-caffeoyl": "*[C@H]1O[C@@H](CO)[C@H](O)[C@@H](OC(=O)/C=C/c2ccc(O)c(O)c2)[C@H]1O",
}


DROP_PLACEHOLDER_IDS = {
    "LIT_0001",
    "LIT_0002",
    "LIT_0126",
    "LIT_0127",
    "LIT_0145",
    "LIT_0147",
    "LIT_0148",
    "LIT_0149",
    "LIT_0150",
    "LIT_0151",
    "LIT_0152",
    "LIT_0153",
    "LIT_0154",
    "LIT_0155",
    "LIT_0261",
    "LIT_0262",
    "LIT_0263",
    "LIT_0264",
    "LIT_0265",
    "LIT_0266",
    "LIT_0267",
    "LIT_0268",
    "LIT_0269",
    "LIT_0281",
    "LIT_0109",
    "LIT_0110",
    "LIT_0111",
    "LIT_0112",
    "LIT_0113",
    "LIT_0114",
    "LIT_0117",
    "LIT_0120",
    "LIT_0124",
    "LIT_0179",
    "LIT_0186",
    "LIT_0202",
    "LIT_0220",
    "LIT_0225",
    "LIT_0226",
    "LIT_0227",
    "LIT_0228",
    "LIT_0229",
    "LIT_0230",
    "LIT_0246",
    "LIT_0247",
    "LIT_0248",
    "LIT_0249",
    "LIT_0334",
}


def normalize_manual_placeholder_text(text: str) -> str:
    return (
        text.replace("[R]", "*")
        .replace("[Ar]", "*")
        .replace("[X]", "*")
        .replace("[Y]", "*")
        .replace("[*:1]", "*")
        .replace("[*:2]", "*")
        .replace("[*]", "*")
    )


def count_placeholders_in_component(component: str) -> int:
    mol = Chem.MolFromSmiles(normalize_manual_placeholder_text(component))
    if mol is None:
        return normalize_manual_placeholder_text(component).count("*")
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0)


def count_placeholders(raw: str) -> int:
    return sum(
        count_placeholders_in_component(component)
        for side in raw.split(">>")
        for component in side.split(".")
        if component.strip()
    )


def fragment_smiles(label: str) -> str:
    label = label.strip()
    if label not in ATTACHMENT_FRAGMENTS:
        raise KeyError(f"No attachment fragment for {label!r}")
    return ATTACHMENT_FRAGMENTS[label]


def substitute_dummy_once(mol: Chem.Mol, dummy_idx: int, fragment: str) -> Chem.Mol:
    dummy = mol.GetAtomWithIdx(dummy_idx)
    neighbors = [atom.GetIdx() for atom in dummy.GetNeighbors()]
    if len(neighbors) != 1:
        raise ValueError("placeholder must have exactly one neighbor")
    base_neighbor = neighbors[0]
    if fragment == "[H]":
        rw = Chem.RWMol(mol)
        rw.RemoveAtom(dummy_idx)
        out = rw.GetMol()
        Chem.SanitizeMol(out)
        return out

    frag = Chem.MolFromSmiles(fragment)
    if frag is None:
        raise ValueError(f"invalid fragment SMILES: {fragment}")
    frag_dummy_indices = [atom.GetIdx() for atom in frag.GetAtoms() if atom.GetAtomicNum() == 0]
    if len(frag_dummy_indices) != 1:
        raise ValueError(f"fragment must have exactly one attachment dummy: {fragment}")
    frag_dummy_idx = frag_dummy_indices[0]
    frag_dummy = frag.GetAtomWithIdx(frag_dummy_idx)
    frag_neighbors = [atom.GetIdx() for atom in frag_dummy.GetNeighbors()]
    if len(frag_neighbors) != 1:
        raise ValueError(f"fragment attachment dummy must have one neighbor: {fragment}")
    frag_neighbor = frag_neighbors[0]

    offset = mol.GetNumAtoms()
    rw = Chem.RWMol(Chem.CombineMols(mol, frag))
    rw.AddBond(base_neighbor, offset + frag_neighbor, Chem.BondType.SINGLE)
    for idx in sorted([dummy_idx, offset + frag_dummy_idx], reverse=True):
        rw.RemoveAtom(idx)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


def substitute_component_placeholders(component: str, fragments: list[str], start: int) -> tuple[str, int]:
    normalized = normalize_manual_placeholder_text(component.strip())
    if "*" not in normalized:
        return component.strip(), start
    mol = Chem.MolFromSmiles(normalized)
    if mol is None:
        raise ValueError(f"cannot parse placeholder component: {component}")
    for fragment in fragments[start : start + count_placeholders_in_component(component)]:
        dummy_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
        if not dummy_indices:
            raise ValueError(f"too many fragments for component: {component}")
        mol = substitute_dummy_once(mol, dummy_indices[0], fragment)
    used = count_placeholders_in_component(component)
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        raise ValueError(f"unsubstituted placeholders remain in component: {component}")
    return Chem.MolToSmiles(mol, isomericSmiles=True), start + used


def substitute_reaction_placeholders(raw: str, labels: list[str]) -> str:
    fragments = [fragment_smiles(label) for label in labels]
    expected = count_placeholders(raw)
    if len(fragments) != expected:
        raise ValueError(f"expected {expected} fragments for {raw}, got {len(fragments)}")
    used = 0
    sides = []
    for side in raw.split(">>"):
        components = []
        for component in [part.strip() for part in side.split(".") if part.strip()]:
            substituted, used = substitute_component_placeholders(component, fragments, used)
            components.append(substituted)
        sides.append(".".join(components))
    if used != len(fragments):
        raise ValueError(f"unused fragments: used {used}, got {len(fragments)}")
    return ">>".join(sides)


def split_meaning_values(text: str) -> list[str]:
    text = text.strip()
    text = re.sub(r"^[：:=\s]+", "", text)
    text = re.sub(r"[。；;，,、\s]+$", "", text)
    text = text.replace("或", ",").replace("和", ",").replace("及", ",")
    text = text.replace("，", ",").replace("、", ",").replace("；", ",").replace(";", ",")
    values = [v.strip(" -()（）`") for v in text.split(",") if v.strip(" -()（）`")]
    normalized = []
    for value in values:
        value = value.replace("-OMe", "OMe").replace("-CH=CH-CH=CH2", "butadienyl")
        value = value.replace("烷基", "alkyl").replace("芳基", "aryl").replace("芳香环", "aryl")
        value = value.replace("卤素", "halogen").replace("葡萄糖", "Glu")
        if value in ("MeO", "OMe"):
            value = "OMe"
        if value in ("β-Glc", "β-Glu", "Gle"):
            value = "Glu"
        normalized.append(value)
    return normalized


def possible_structures(labels: list[str]) -> list[dict[str, Any]]:
    structures = []
    for label in labels:
        key = label
        if key.startswith("Glu-") or key.startswith("Glc-"):
            base = dict(SUBSTITUENT_LIBRARY["Glu"])
            base["label"] = label
            base["type"] = "decorated_" + base["type"]
            base["decoration"] = label
            structures.append(base)
            continue
        if key.startswith("Rha-"):
            base = dict(SUBSTITUENT_LIBRARY["Rha"])
            base["label"] = label
            base["type"] = "decorated_" + base["type"]
            base["decoration"] = label
            structures.append(base)
            continue
        if key in SUBSTITUENT_LIBRARY:
            base = dict(SUBSTITUENT_LIBRARY[key])
            base["label"] = label
            structures.append(base)
        elif re.fullmatch(r"\*?C{3,10}|C{3,10}", key):
            structures.append({"label": label, "type": "alkyl_chain_or_unsaturated_chain", "attachment_smarts": "[*:1][C]", "smiles": ""})
        elif key in ("butadienyl",):
            structures.append({"label": label, "type": "dienyl", "attachment_smarts": "[*:1]/C=C/C=C", "smiles": "[*:1]/C=C/C=C"})
        else:
            structures.append({"label": label, "type": "unmapped_textual_meaning", "attachment_smarts": "", "smiles": ""})
    return structures


def find_symbol_meanings(raw: str, context: str) -> list[dict[str, Any]]:
    annotations: list[dict[str, Any]] = []
    symbols = []
    for symbol in ("R1", "R2", "R", "Ar", "X", "Y", "*"):
        if symbol == "*":
            if "*" in raw or "\\*" in context or "* =" in context:
                symbols.append(symbol)
        elif f"[{symbol}]" in raw or re.search(rf"\b{re.escape(symbol)}\s*[=：:]", context):
            symbols.append(symbol)

    def add(symbol: str, labels: list[str], evidence: str, site: str = "", exactness: str = "finite") -> None:
        labels = [label for label in labels if label]
        if not labels:
            return
        annotations.append(
            {
                "symbol": symbol,
                "site": site,
                "meanings": labels,
                "possible_structures": possible_structures(labels),
                "evidence": evidence.strip(),
                "exactness": exactness,
            }
        )

    # Site-specific sugar annotations in the notes.
    for site, value in re.findall(r"([ABC]环)\\\\?\*\\s*=\\s*([^，；;。`\\n：]+)", context):
        add("*", split_meaning_values(value), f"{site}* = {value}", site=site, exactness="exact")
    for phase, value in re.findall(r"(反应物|产物)的\\*\\s*=\\s*([^，；;。`\\n]+)", context):
        add("*", split_meaning_values(value), f"{phase}* = {value}", site=phase, exactness="exact")
    for site, value in re.findall(r"([ABC]环)?\\s*\\*\\s*=\\s*([^；;。`\\n]+)", context):
        add("*", split_meaning_values(value), f"{site or ''}* = {value}", site=site or "", exactness="exact")
    for site, value in re.findall(r"([ABC]环)\\s*\\*\\s*=\\s*([^；;。`\\n]+)", context):
        add("*", split_meaning_values(value), f"{site}* = {value}", site=site, exactness="exact")
    for site, value in re.findall(r"([ABC]环)\\\\\\*\\s*=\\s*([^；;。`\\n]+)", context):
        add("*", split_meaning_values(value), f"{site}* = {value}", site=site, exactness="exact")

    # Generic assignments such as R = -SCoA, R1 = OH, X = H, Cl...
    for symbol, value in re.findall(r"\\b(R1|R2|R|X|Y|Ar)\\s*[=：:]\\s*([^。；;`\\n]+)", context):
        add(symbol, split_meaning_values(value), f"{symbol} = {value}", exactness="finite")

    if "R = -SCoA" in context or "R = - SCoA" in context:
        add("R", ["SCoA"], "R = -SCoA", exactness="exact")
    if "Ar表示芳香环" in context or "Ar：芳基" in context or "Ar, R：芳基" in context:
        add("Ar", ["aryl"], "Ar is defined as an aromatic/aryl group", exactness="broad_class")
    if "Ar, R：芳基" in context or "Ar和R均为芳基" in context:
        add("R", ["aryl"], "R is defined as aryl in this context", exactness="broad_class")
    if "R 为 H 或烷基" in context or "R为 H 或烷基" in context:
        add("R", ["H", "alkyl"], "R is H or alkyl", exactness="broad_class")
    if "X = H, Cl, Br, Me, OMe" in context:
        add("X", ["H", "Cl", "Br", "Me", "OMe"], "X = H, Cl, Br, Me, OMe", exactness="finite")
    if "Y = H, Me, OMe" in context:
        add("Y", ["H", "Me", "OMe"], "Y = H, Me, OMe", exactness="finite")
    if "R = H，烷基，烷氧基" in context or "R = H, 烷基, 烷氧基" in context:
        add("R", ["H", "alkyl", "OMe"], "R = H, alkyl, alkoxy", exactness="broad_class")
    if "* = Enz" in context or "* =  Enz" in context:
        add("*", ["Enz"], "* = Enz", exactness="cofactor")
    if "* = OPP" in context or "(* = OPP)" in context:
        add("*", ["OPP"], "* = OPP", exactness="exact")
    if "* = Glu" in context or "下面的* = Glu" in context:
        add("*", ["Glu"], "* = Glu", exactness="exact")
    if "A环\\* = Rha" in context:
        add("*", ["Rha"], "A-ring * = Rha", site="A环", exactness="exact")
    if "C环\\* = Api" in context:
        add("*", ["Api"], "C-ring * = Api", site="C环", exactness="exact")
    if "C环\\* = Glc(3-Glc)-4-caffeoyl" in context:
        add("*", ["Glc(3-Glc)-4-caffeoyl"], "C-ring * = Glc(3-Glc)-4-caffeoyl", site="C环", exactness="exact")
    if "反应物的* = Glc(3-Glc)-4-caffeoyl" in context:
        add("*", ["Glc(3-Glc)-4-caffeoyl"], "reactant * = Glc(3-Glc)-4-caffeoyl", site="reactant", exactness="exact")
    if "产物的* = Glc-caffeoyl" in context:
        add("*", ["Glc-caffeoyl"], "product * = Glc-caffeoyl", site="product", exactness="exact")
    if "Glu-(2-1)-Rha" in context:
        add("*", ["Glu-(2-1)-Rha"], "* = Glu-(2-1)-Rha", exactness="exact")
    if "R1 = OH，R2 = H" in context:
        add("R1/R2", ["R1:OH;R2:H", "R1:OMe;R2:OH"], "Naringin/neohesperidin finite variants", exactness="finite")
    if "* = H, -OMe" in context or "* = H, OMe" in context:
        add("*", ["H", "OMe"], "* = H, OMe", exactness="finite")
    if "* = H, MeO" in context:
        add("*", ["H", "OMe"], "* = H, MeO", exactness="finite")
    if "R = H, NO2, MeO, CF3" in context or "R = H，NO2，MeO，CF3" in context:
        add("R", ["H", "NO2", "OMe", "CF3"], "R = H, NO2, MeO, CF3", exactness="finite")
    if "R = `$smiles=*[Si](C)(C)C`" in context:
        add("R", ["trimethylsilyl", "hexyl", "octyl", "allyl", "butenyl", "styryl"], "R examples are listed as TMS, hexyl, octyl, allyl, butenyl, styryl SMILES", exactness="finite")
    if "* = H，F" in context:
        add("*", ["H", "F"], "* = H, F", exactness="finite")
    if "(R = H, F)" in context:
        add("R", ["H", "F"], "R = H, F", exactness="finite")
    if "* = H, Me; R = H, OMe" in context:
        add("*", ["H", "Me"], "* = H, Me", exactness="finite")
        add("R", ["H", "OMe", "i-Pr"], "R = H, OMe, i-Pr", exactness="finite")
    if "* = H, OMe" in context:
        add("*", ["H", "OMe"], "* = H, OMe", exactness="finite")

    # If a placeholder remains but context only says substituted ring, record that
    # explicitly instead of silently treating it as any atom.
    if not annotations and symbols:
        aryl_placeholder = any(token in raw for token in ("[*]c", "c([*])", "c(*)", "(*)", "[R]c", "c([R])"))
        if "任意位点取代基" in context or "取代基" in context or aryl_placeholder:
            for symbol in symbols:
                add(symbol, ["aryl_ring_substituent"], "Context only states substituted aryl/ring substituent", exactness="broad_class")
        else:
            for symbol in symbols:
                add(symbol, ["unresolved_placeholder"], "No explicit meaning found in nearby context", exactness="unresolved")

    # Deduplicate annotations.
    deduped = []
    seen = set()
    for ann in annotations:
        key = (ann["symbol"], ann["site"], tuple(ann["meanings"]), ann["evidence"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(ann)
    return deduped


def generic_resolution_status(annotations: list[dict[str, Any]], has_generic: bool) -> str:
    if not has_generic:
        return "not_generic"
    if not annotations:
        return "unresolved"
    exactness = {ann.get("exactness", "") for ann in annotations}
    if "unresolved" in exactness:
        return "unresolved"
    if "cofactor" in exactness:
        return "cofactor_or_enzyme_placeholder"
    if "broad_class" in exactness:
        return "broad_class"
    if exactness <= {"exact", "finite"}:
        return "resolved_finite_or_exact"
    return "partially_resolved"


@dataclass
class RawReaction:
    raw: str
    source_file: Path
    start_line: int
    extraction_method: str
    heading_path: list[str]
    context_before: list[str]
    context_after: list[str]


def canonical_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def inchikey(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return ""


def normalize_placeholder_text(text: str) -> str:
    replacements = {
        "[R]": "[*]",
        "[r]": "[*]",
        "[Ar]": "[*]",
        "[AR]": "[*]",
        "[X]": "[F,Cl,Br,I]",
        "[Y]": "[F,Cl,Br,I]",
        "[HH]": "[H][H]",
        "O=[Fe]-*": "O=[Fe]-[*]",
        "O[Fe]*": "O[Fe][*]",
        "[Fe]=*": "[Fe]=[*]",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text.strip()


def component_to_smarts(component: str) -> tuple[str, bool, str]:
    component = normalize_placeholder_text(component.strip())
    if not component:
        return "", False, "empty component"
    mol = Chem.MolFromSmiles(component)
    if mol is not None:
        return Chem.MolToSmarts(mol, isomericSmiles=True), True, ""
    patt = Chem.MolFromSmarts(component)
    if patt is not None:
        return component, True, "kept as SMARTS"
    return component, False, "invalid SMILES/SMARTS component"


def side_to_smarts(side: str) -> tuple[str, bool, list[str]]:
    parts = [p.strip() for p in side.split(".") if p.strip()]
    converted: list[str] = []
    errors: list[str] = []
    ok = True
    for part in parts:
        smarts, part_ok, note = component_to_smarts(part)
        converted.append(smarts)
        if not part_ok:
            ok = False
            errors.append(f"{part}: {note}")
    return ".".join(converted), ok, errors


def validate_reaction_smarts(smarts: str) -> tuple[bool, str]:
    if ">>" not in smarts:
        return False, "no reaction arrow"
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        if rxn is None:
            return False, "RDKit returned None"
        if rxn.GetNumReactantTemplates() == 0 or rxn.GetNumProductTemplates() == 0:
            return False, (
                f"reactants={rxn.GetNumReactantTemplates()}, "
                f"products={rxn.GetNumProductTemplates()}"
            )
        return True, (
            f"reactants={rxn.GetNumReactantTemplates()}, "
            f"products={rxn.GetNumProductTemplates()}"
        )
    except Exception as exc:
        return False, str(exc)[:220]


def reaction_to_smarts(raw: str) -> dict[str, Any]:
    left, right = [part.strip() for part in raw.split(">>", 1)]
    left_smarts, left_ok, left_errors = side_to_smarts(left)
    right_smarts, right_ok, right_errors = side_to_smarts(right)
    forward = f"{left_smarts}>>{right_smarts}"
    retro = f"{right_smarts}>>{left_smarts}"
    forward_ok, forward_msg = validate_reaction_smarts(forward)
    retro_ok, retro_msg = validate_reaction_smarts(retro)
    return {
        "reactants_smiles": left,
        "products_smiles": right,
        "forward_smarts": forward,
        "retro_template": retro,
        "component_conversion_ok": left_ok and right_ok,
        "component_errors": left_errors + right_errors,
        "forward_validation": {"ok": forward_ok, "message": forward_msg},
        "retro_validation": {"ok": retro_ok, "message": retro_msg},
    }


def split_components(side: str) -> list[dict[str, Any]]:
    output = []
    for component in [p.strip() for p in side.split(".") if p.strip()]:
        cleaned = normalize_placeholder_text(component)
        has_dummy = "*" in cleaned or "[*" in cleaned or "[R]" in component or "[Ar]" in component
        can = ""
        key = ""
        mw = None
        if not has_dummy:
            can = canonical_smiles(cleaned)
            key = inchikey(can) if can else ""
            if can:
                mol = Chem.MolFromSmiles(can)
                mw = round(Descriptors.MolWt(mol), 3) if mol is not None else None
        output.append(
            {
                "raw": component,
                "normalized": cleaned,
                "canonical_smiles": can,
                "inchikey": key,
                "molecular_weight": mw,
                "has_wildcard_or_generic_group": has_dummy,
            }
        )
    return output


def update_headings(headings: list[tuple[int, str]], line: str) -> list[tuple[int, str]]:
    match = re.match(r"^(#{1,6})\s+(.+?)\s*$", line)
    if not match:
        return headings
    level = len(match.group(1))
    title = match.group(2).strip()
    headings = [(old_level, old_title) for old_level, old_title in headings if old_level < level]
    headings.append((level, title))
    return headings


def clean_reaction_candidate(raw: str) -> str:
    raw = raw.strip().strip("`").strip()
    raw = raw.replace("&gt;&gt;", ">>")
    raw = re.sub(r"^\$smiles=", "", raw)
    raw = re.sub(r"[，。,;；：:、]+$", "", raw)
    # One source reaction has a truncated p-iodophenyl acetate reactant.
    raw = raw.replace("[I]c1ccc(OC(=O)C)>>", "[I]c1ccc(OC(=O)C)cc1>>")
    return raw


def extract_reactions_from_file(path: Path) -> list[RawReaction]:
    lines = path.read_text(encoding="utf-8").splitlines()
    headings: list[tuple[int, str]] = []
    reactions: list[RawReaction] = []
    in_smiles_block = False
    block_start = 0
    block_headings: list[tuple[int, str]] = []
    block_lines: list[tuple[int, str]] = []

    def context(line_no: int) -> tuple[list[str], list[str]]:
        start = max(0, line_no - 9)
        end = min(len(lines), line_no + 6)
        before = [line.strip() for line in lines[start : max(0, line_no - 1)] if line.strip()]
        after = [line.strip() for line in lines[max(0, line_no - 1) : end] if line.strip()]
        return before, after

    for idx, line in enumerate(lines, start=1):
        headings = update_headings(headings, line)
        stripped = line.strip()

        if stripped.startswith("```smiles"):
            in_smiles_block = True
            block_start = idx
            block_headings = list(headings)
            block_lines = []
            continue

        if in_smiles_block and stripped.startswith("```"):
            for raw_line_no, raw_line in block_lines:
                candidate = clean_reaction_candidate(raw_line)
                if ">>" not in candidate:
                    continue
                before, after = context(raw_line_no)
                reactions.append(
                    RawReaction(
                        raw=candidate,
                        source_file=path,
                        start_line=raw_line_no,
                        extraction_method=f"smiles_code_block_start_{block_start}",
                        heading_path=[title for _, title in block_headings],
                        context_before=before,
                        context_after=after,
                    )
                )
            in_smiles_block = False
            continue

        if in_smiles_block:
            block_lines.append((idx, stripped))
            continue

        for match in re.finditer(r"\$smiles=([^`\s]+)", line):
            candidate = clean_reaction_candidate(match.group(1))
            if ">>" not in candidate:
                continue
            before, after = context(idx)
            reactions.append(
                RawReaction(
                    raw=candidate,
                    source_file=path,
                    start_line=idx,
                    extraction_method="inline_smiles",
                    heading_path=[title for _, title in headings],
                    context_before=before,
                    context_after=after,
                )
            )
    return reactions


def parse_conditions(raw: RawReaction) -> dict[str, Any]:
    context = raw.context_before + raw.context_after
    condition_lines = [
        line
        for line in context
        if any(keyword in line for keyword in CONDITION_KEYWORDS)
        or re.search(r"\b(?:Ag2O|Ag2CO3|BF3|SnCl4|TMSOTf|TBAB|DDQ|Pd/C|K2CO3)\b", line)
    ]
    condition_text = " ".join(condition_lines)

    def collect(patterns: dict[str, str]) -> list[str]:
        found = []
        for label, pattern in patterns.items():
            if re.search(pattern, condition_text, flags=re.IGNORECASE):
                found.append(label)
        return sorted(set(found))

    yields = re.findall(r"(?:产率|收率)[^\d]{0,12}(\d+(?:\.\d+)?\s*%)(?:\s*[-~–]\s*(\d+(?:\.\d+)?\s*%))?", condition_text)
    yields_flat = ["-".join([part for part in pair if part]) for pair in yields]
    if not yields_flat:
        yields_flat = re.findall(r"\b\d+(?:\.\d+)?\s*%", condition_text)

    temperatures = re.findall(r"-?\d+(?:\.\d+)?\s*(?:°C|℃)|室温|低温|回流", condition_text)
    time_values = re.findall(r"\b\d+(?:\.\d+)?\s*(?:h|hr|小时|min|分钟|d|天)\b", condition_text, flags=re.IGNORECASE)
    operations = [
        word
        for word in ("回流", "加热", "低温", "室温", "搅拌", "脱保护", "氧化", "还原", "环化", "糖苷化")
        if word in condition_text
    ]

    caption = ""
    for line in raw.context_after:
        if re.match(r"^(图|Scheme|路线)\s*", line):
            caption = line
            break

    return {
        "condition_text": condition_text,
        "condition_lines": condition_lines,
        "catalysts_or_promoters": collect(CATALYST_PATTERNS),
        "solvents": collect(SOLVENT_PATTERNS),
        "temperatures": sorted(set(temperatures)),
        "times": sorted(set(time_values)),
        "yields": sorted(set(yields_flat)),
        "operations": sorted(set(operations)),
        "nearby_caption": caption,
        "context_before": raw.context_before[-5:],
        "context_after": raw.context_after[:5],
    }


def classify_reaction(raw: RawReaction, converted: dict[str, Any], conditions: dict[str, Any]) -> str:
    text = " ".join(raw.heading_path + conditions.get("condition_lines", []) + [raw.raw]).lower()
    chinese_text = " ".join(raw.heading_path + conditions.get("condition_lines", []))
    reactants = converted["reactants_smiles"]
    products = converted["products_smiles"]
    full = f"{reactants}>>{products}"

    if "生物合成" in raw.source_file.name or "生物合成" in chinese_text:
        return "biosynthesis"
    if "糖苷" in chinese_text or "glycos" in text:
        if re.search(r"\bBrC|\bFC|O=C\(NC\(Cl\)\(Cl\)Cl\)", reactants):
            return "o_glycosylation"
        if re.search(r"\[li\]|\[Li\]|傅-克|friedel", full + chinese_text, flags=re.IGNORECASE):
            return "c_glycosylation"
        return "glycosylation"
    if "脱保护" in chinese_text or re.search(r"OCC\d*=CC=CC=C\d*|OCc\d*ccccc", reactants) and ".O" in full:
        return "deprotection"
    if "保护" in chinese_text or "OCOC" in products or "OC(C)=O" in products:
        return "protection"
    if "查耳酮" in chinese_text or "chalcone" in text or "/C=C/" in products and "C(=O)" in products:
        return "chalcone_synthesis"
    if "异黄酮" in chinese_text or "isoflav" in text:
        return "isoflavone_synthesis"
    if "二氢黄酮" in chinese_text or "flavanone" in text:
        return "flavanone_synthesis"
    if "黄酮" in chinese_text or "flavone" in text:
        if "环化" in chinese_text or "cycl" in text:
            return "flavonoid_cyclization"
        return "flavonoid_synthesis"
    if "环化" in chinese_text or "cycl" in text:
        return "cyclization"
    if "氧化" in chinese_text or any(token in full for token in ("[O]", "DDQ", "Mn", "O=O")):
        return "oxidation"
    if "还原" in chinese_text or any(token in full for token in ("[H-]", "[BH4-]", "[AlH4-]", "[H][H]")):
        return "reduction"
    if "甲基化" in chinese_text or "methyl" in text or "CI" in reactants:
        return "methylation"
    if "prenyl" in text or "异戊烯" in chinese_text or "CC=C(C)C" in products:
        return "prenylation"
    return "literature_reaction"


def confidence_for(entry: dict[str, Any]) -> float:
    score = 0.35
    if entry["validation"]["retro"]["ok"]:
        score += 0.25
    if entry["validation"]["component_conversion_ok"]:
        score += 0.15
    if entry["conditions"]["condition_text"]:
        score += 0.10
    if entry["conditions"]["nearby_caption"]:
        score += 0.05
    if not entry["flags"]["has_generic_group"]:
        score += 0.10
    return round(min(score, 0.95), 2)


def rebuild_entry_from_raw(
    source_entry: dict[str, Any],
    raw: str,
    reaction_id: str,
    note: str,
    variant_label: str = "",
) -> dict[str, Any]:
    converted = reaction_to_smarts(raw)
    reactant_components = split_components(converted["reactants_smiles"])
    product_components = split_components(converted["products_smiles"])
    has_generic = any(c["has_wildcard_or_generic_group"] for c in reactant_components + product_components)
    generic_context = " ".join(
        source_entry["conditions"].get("context_before", [])
        + source_entry["conditions"].get("context_after", [])
        + source_entry["conditions"].get("condition_lines", [])
        + [note]
    )
    generic_groups = find_symbol_meanings(raw, generic_context) if has_generic else []
    resolution_status = generic_resolution_status(generic_groups, has_generic)

    entry = copy.deepcopy(source_entry)
    entry["reaction_id"] = reaction_id
    entry["raw_reaction_smiles"] = raw
    entry["reactants"] = reactant_components
    entry["products"] = product_components
    entry["reaction_smarts"] = {
        "forward": converted["forward_smarts"],
        "retro": converted["retro_template"],
    }
    entry["generic_groups"] = generic_groups
    entry["validation"] = {
        "component_conversion_ok": converted["component_conversion_ok"],
        "component_errors": converted["component_errors"],
        "forward": converted["forward_validation"],
        "retro": converted["retro_validation"],
    }
    entry["flags"] = {
        **entry["flags"],
        "has_generic_group": has_generic,
        "generic_resolution_status": resolution_status,
        "aizynthfinder_ready": converted["retro_validation"]["ok"] and not has_generic,
        "has_multiple_products": len(product_components) > 1,
        "has_multiple_reactants": len(reactant_components) > 1,
    }
    entry["curation"] = {
        **entry["curation"],
        "manual_review_applied": True,
        "parent_reaction_id": source_entry["reaction_id"],
        "variant_label": variant_label,
        "manual_review_note": note,
        "needs_review": bool(has_generic or not converted["retro_validation"]["ok"]),
        "review_reasons": [],
    }
    if has_generic:
        entry["curation"]["review_reasons"].append(f"contains wildcard/generic group ({resolution_status})")
    if not converted["retro_validation"]["ok"]:
        entry["curation"]["review_reasons"].append(converted["retro_validation"]["message"])
    if note:
        entry["conditions"]["manual_review_note"] = note
        if "enzyme" in note.lower() or "酶" in note:
            entry["conditions"]["catalysts_or_promoters"] = sorted(
                set(entry["conditions"].get("catalysts_or_promoters", []) + ["enzyme/cofactor"])
            )
    entry["curation"]["confidence"] = confidence_for(entry)
    return entry


def manual_variant(
    entry: dict[str, Any],
    raw: str | None,
    labels: list[str] | None,
    note: str,
    label: str = "",
) -> dict[str, Any]:
    source_raw = raw or entry["raw_reaction_smiles"]
    curated_raw = substitute_reaction_placeholders(source_raw, labels) if labels is not None else source_raw
    return rebuild_entry_from_raw(entry, curated_raw, entry["reaction_id"], note, label)


def same_value_variants(entry: dict[str, Any], values: list[str], note: str) -> list[dict[str, Any]]:
    raw = entry["raw_reaction_smiles"]
    n_placeholders = count_placeholders(raw)
    variants = []
    for value in values:
        variants.append(
            manual_variant(
                entry,
                raw,
                [value] * n_placeholders,
                note,
                value,
            )
        )
    return variants


def with_variant_ids(entry: dict[str, Any], variants: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if len(variants) == 1:
        variants[0]["reaction_id"] = entry["reaction_id"]
        return variants
    for idx, variant in enumerate(variants, start=1):
        variant["reaction_id"] = f"{entry['reaction_id']}_v{idx:02d}"
    return variants


def replace_component(raw: str, old: str, new: str) -> str:
    if old not in raw:
        raise ValueError(f"component not found: {old}")
    return raw.replace(old, new, 1)


def manual_curated_variants(entry: dict[str, Any]) -> list[dict[str, Any]] | None:
    rid = entry["reaction_id"]
    raw = entry["raw_reaction_smiles"]

    if rid == "LIT_0109":
        fixed = "c1ccc(cc1)C(=O)C[Te+](CCCC)(CCCC).[Br-].O=Cc1ccccc1>>c1ccc(cc1)C(=O)C=Cc1ccccc1"
        return with_variant_ids(
            entry,
            [
                rebuild_entry_from_raw(
                    entry,
                    fixed,
                    entry["reaction_id"],
                    "Manual review: broad R/Ar collapsed to the unsubstituted aryl skeleton; Te/Br ion pair disconnected for RDKit-valid template export.",
                    "skeleton",
                )
            ],
        )

    if rid == "LIT_0117":
        fixed = "Clc1ccccc1.C#CC(O)c1ccccc1>>c1ccccc1C=CC(O)c1ccccc1"
        return with_variant_ids(
            entry,
            [
                rebuild_entry_from_raw(
                    entry,
                    fixed,
                    entry["reaction_id"],
                    "Manual review: broad Ar/X collapsed to the parent chlorophenyl skeleton.",
                    "skeleton",
                )
            ],
        )

    if rid in DROP_PLACEHOLDER_IDS:
        return with_variant_ids(entry, same_value_variants(entry, ["H"], "Manual review: placeholders were removed/replaced by H."))

    if rid in {"LIT_0012", "LIT_0013", "LIT_0015", "LIT_0017", "LIT_0018"}:
        return with_variant_ids(
            entry,
            [
                manual_variant(
                    entry,
                    raw,
                    ["SCoA"] * count_placeholders(raw),
                    "Manual review: R = SCoA; encoded as a thioester sulfur attachment for template use.",
                    "SCoA",
                )
            ],
        )

    if rid in {"LIT_0020", "LIT_0021"}:
        variants = []
        for variant_name, r2, r1 in (
            ("R1=OH_R2=H", "H", "OH"),
            ("R1=OMe_R2=OH", "OH", "OMe"),
        ):
            variants.append(
                manual_variant(
                    entry,
                    raw,
                    ["Glc-Rha", r2, r1, "Glc-Rha", r2, r1],
                    "Manual review: * = Glc-Rha; first [R] is R2 and second [R] is R1.",
                    variant_name,
                )
            )
        return with_variant_ids(entry, variants)

    if rid == "LIT_0022":
        fixed = "O=C1C[C@H](c2ccc(O)cc2)Oc2cc(O)cc([R])c21>>O=C1[CH][C@H](c2ccc(O)cc2)Oc2cc(O)cc([R])c21"
        variants = [
            manual_variant(entry, fixed, [value, value], "Manual review: Fe-Enz cofactor removed from SMILES and recorded as enzyme/cofactor condition.", value)
            for value in ("H", "OH")
        ]
        return with_variant_ids(entry, variants)

    if rid == "LIT_0024":
        fixed = "O=C1[CH][C@H](c2ccc(O)cc2)Oc2cc(O)cc([R])c21>>O=C1[C@@H](c2ccc(O)cc2)[CH]Oc2cc(O)cc([R])c21"
        variants = [
            manual_variant(entry, fixed, [value, value], "Manual review: corrected 1,2-aryl migration; enzyme/cofactor kept in conditions.", value)
            for value in ("H", "OH")
        ]
        return with_variant_ids(entry, variants)

    if rid == "LIT_0025":
        fixed = "O=C1[C@@H](c2ccc(O)cc2)[CH]Oc2cc(O)cc([R])c21>>O=C1[C@@H](c2ccc(O)cc2)[C@@H](O)Oc2cc(O)cc([R])c21"
        variants = [
            manual_variant(entry, fixed, [value, value], "Manual review: Fe-Enz cofactor removed from SMILES and recorded as enzyme/cofactor condition.", value)
            for value in ("H", "OH")
        ]
        return with_variant_ids(entry, variants)

    if rid == "LIT_0034":
        fixed = raw.replace("CC(=C)=CC*", "CC(C)=CCOP(=O)([O-])OP(=O)([O-])[O-]")
        return with_variant_ids(entry, [manual_variant(entry, fixed, None, "Manual review: * = OPP; encoded as dimethylallyl diphosphate.", "OPP")])

    if rid == "LIT_0088":
        return with_variant_ids(entry, [manual_variant(entry, raw, ["Rha", "Api", "Rha"], "Manual review: first reactant * = Rha, second reactant * = Api, product * = Rha.", "Rha_Api")])
    if rid == "LIT_0089":
        return with_variant_ids(entry, [manual_variant(entry, raw, ["Rha"], "Manual review: reactant * = Rha.", "Rha")])
    if rid == "LIT_0090":
        return with_variant_ids(entry, [manual_variant(entry, raw, ["Rha", "Glc(3-Glc)-4-caffeoyl", "Glc(3-Glc)-4-caffeoyl"], "Manual review: reactant *=Rha and Glc(3-Glc)-4-caffeoyl; product *=Glc(3-Glc)-4-caffeoyl.", "Rha_GlcGlcCaff")])
    if rid == "LIT_0091":
        return with_variant_ids(entry, [manual_variant(entry, raw, ["Glc(3-Glc)-4-caffeoyl", "Glc-caffeoyl"], "Manual review: reactant * = Glc(3-Glc)-4-caffeoyl; product * = Glc-caffeoyl.", "GlcGlcCaff_to_GlcCaff")])
    if rid == "LIT_0094":
        return with_variant_ids(entry, [manual_variant(entry, raw, ["Glu", "Glu"], "Manual review: * = Glu.", "Glu")])

    if rid == "LIT_0101":
        variants = []
        for first, halogen, second in itertools.product(("H", "Me", "F", "CF3"), ("Cl", "Br", "I"), ("CF3", "NO2", "Me")):
            variants.append(manual_variant(entry, raw, [first, halogen, second], "Manual review: finite substituent expansion for *, X.", f"{first}_{halogen}_{second}"))
        return with_variant_ids(entry, variants)

    if rid == "LIT_0105":
        variants = []
        for star, r_value in itertools.product(("H", "OMe"), ("H", "NO2", "OMe", "CF3")):
            variants.append(manual_variant(entry, raw, [star, r_value, r_value, star], "Manual review: finite expansion for * and R.", f"star={star}_R={r_value}"))
        return with_variant_ids(entry, variants)

    if rid == "LIT_0106":
        variants = []
        for r_value, halogen in itertools.product(("TMS", "hexyl", "octyl", "allyl", "butenyl", "styryl"), ("Cl", "Br", "I")):
            variants.append(manual_variant(entry, raw, [r_value, halogen, r_value], "Manual review: X = Cl/Br/I; R expanded from listed substituent SMILES.", f"{r_value}_{halogen}"))
        return with_variant_ids(entry, variants)

    if rid == "LIT_0107":
        x_scaffolds = [
            ("H", "COC(=O)C=Cc1ccccc1", "C=Cc1ccccc1"),
            ("2-OMe", "COc1ccccc1C=CC(=O)OC", "c1ccccc1OC"),
            ("2,4-(OMe)2", "COc1ccc(OC)cc1C=CC(=O)OC", "c1cc(OC)ccc1OC"),
            ("2,6-(OMe)2", "COc1cccc(OC)c1C=CC(=O)OC", "c1cccc(OC)c1OC"),
            ("2,4,6-(OMe)3", "COc1cc(OC)cc(OC)c1C=CC(=O)OC", "c1cc(OC)cc(OC)c1OC"),
        ]
        y_scaffolds = [
            ("H", "[Li]c1ccccc1", "c1ccccc1"),
            ("4-OMe", "[Li]c1ccc(OC)cc1", "c1ccc(OC)cc1"),
        ]
        variants = []
        for x_label, ester, x_aryl_after_alkene in x_scaffolds:
            for y_label, aryl_li, y_aryl in y_scaffolds:
                if x_label == "H":
                    product = f"O=C({x_aryl_after_alkene}){y_aryl}"
                else:
                    product = ester.replace("C(=O)OC", f"C(=O){y_aryl}")
                fixed = f"{ester}.{aryl_li}>>{product}"
                variants.append(
                    rebuild_entry_from_raw(
                        entry,
                        fixed,
                        entry["reaction_id"],
                        "Manual review: X/Y finite expansion; broad ester R=alkyl represented as methyl for a concrete template.",
                        f"X={x_label}_Y={y_label}_methyl_ester",
                    )
                )
        tin_raw = "c1cc([X])ccc1C=CC(=O)[Cl].[Sn](CCCC)(CCCC)(CCCC)c1ccc([Y])cc1>>c1cc([X])ccc1C=CC(=O)c2ccc([Y])cc2"
        for x_value, y_value in itertools.product(("H", "Cl", "Br", "Me", "OMe"), ("H", "Me", "OMe")):
            variants.append(
                manual_variant(
                    entry,
                    tin_raw,
                    [x_value, y_value, x_value, y_value],
                    "Manual review: added missing organotin coupling; conditions NiCl2(PPh3)3.",
                    f"tin_X={x_value}_Y={y_value}",
                )
            )
        return with_variant_ids(entry, variants)

    if rid in {"LIT_0121", "LIT_0122"}:
        variants = [
            manual_variant(entry, raw, [value] * count_placeholders(raw), "Manual review: R = H or OMe.", value)
            for value in ("H", "OMe")
        ]
        return with_variant_ids(entry, variants)

    if rid == "LIT_0123":
        variants = []
        for r_value, star in itertools.product(("H", "OMe"), ("H", "F")):
            variants.append(manual_variant(entry, raw, [r_value, r_value, star, star, r_value, star, star, r_value], "Manual review: R = H/OMe and * = H/F.", f"R={r_value}_star={star}"))
        return with_variant_ids(entry, variants)

    if rid == "LIT_0129":
        return with_variant_ids(entry, same_value_variants(entry, ["H", "OH"], "Manual review: * = H or OH."))

    if rid == "LIT_0125":
        return with_variant_ids(entry, same_value_variants(entry, ["H", "OMe"], "Manual review: * = H or OMe."))

    if rid == "LIT_0131":
        ketones = [
            "O=C(C#C)c2ccccc2",
            "O=C(C#C)c2ccc(OC)cc2",
            "O=C(C#C)c2c(OC)c(OC)c(OC)cc2",
            "O=C(C#C)c2ccc(O)cc2",
            "O=C(C#C)c2cc(O)c(O)cc2",
            "O=C(C#C)c2cc(O)c(O)c(O)c2",
        ]
        variants = []
        for star, ketone in itertools.product(("H", "Me"), ketones):
            fixed = replace_component(raw, "O=C(C#C)c2ccccc2", ketone)
            variants.append(manual_variant(entry, fixed, [star, star, star, star], "Manual review: * = H/Me and aryl ethynyl ketone expanded from listed reactants.", f"{star}_{ketone}"))
        return with_variant_ids(entry, variants)

    if rid == "LIT_0156":
        variants = same_value_variants(entry, ["H", "OH", "OMe", "Cl", "NO2"], "Manual review: * = H/OH/OMe/Cl/NO2.")
        variants.append(
            rebuild_entry_from_raw(
                entry,
                "Oc1c(C(CC(c2ccccc2)=O)=O)cccc1>>O=c1c2ccccc2oc(c2ccccc2)c1",
                f"{rid}_added_mw",
                "Manual review: added missing microwave cyclization reaction from user note.",
                "added_mw",
            )
        )
        return with_variant_ids(entry, variants)

    if rid in {"LIT_0157", "LIT_0158"}:
        variants = []
        for star, r_value in (("H", "OMe"), ("H", "Me"), ("Me", "Me")):
            variants.append(manual_variant(entry, raw, [star, r_value, star, r_value], "Manual review: three paired * and R combinations.", f"star={star}_R={r_value}"))
        return with_variant_ids(entry, variants)

    if rid in {"LIT_0167", "LIT_0168"}:
        variants = []
        for r_value, star in itertools.product(("H", "Me"), ("H", "Me", "OMe", "Cl")):
            variants.append(manual_variant(entry, raw, [r_value, star, r_value, star], "Manual review: * = H/Me/OMe/Cl; R = H/Me.", f"R={r_value}_star={star}"))
        return with_variant_ids(entry, variants)

    if rid in {"LIT_0231", "LIT_0232"}:
        variants = []
        if rid == "LIT_0231":
            for star, r_value in itertools.product(("H", "OMe"), ("H", "butadienyl")):
                variants.append(manual_variant(entry, raw, [star, star, r_value, r_value], "Manual review: * = H/OMe; R = H/-CH=CH-CH=CH2.", f"star={star}_R={r_value}"))
        else:
            for star, r_value in itertools.product(("H", "OMe"), ("H", "butadienyl")):
                variants.append(manual_variant(entry, raw, [star, r_value, r_value, star, r_value, r_value], "Manual review: * = H/OMe; R = H/-CH=CH-CH=CH2.", f"star={star}_R={r_value}"))
        return with_variant_ids(entry, variants)

    if rid in {"LIT_0256", "LIT_0257", "LIT_0258", "LIT_0259"}:
        return with_variant_ids(entry, same_value_variants(entry, ["H", "Me", "Cl"], "Manual review: all * = H/Me/Cl."))

    if rid in {"LIT_0297", "LIT_0298", "LIT_0299"}:
        return with_variant_ids(entry, same_value_variants(entry, ["H", "F"], "Manual review: R = H/F."))

    if rid in {"LIT_0300", "LIT_0301", "LIT_0302"}:
        variants = []
        for r_value, star in itertools.product(("H", "OMe", "i-Pr"), ("H", "Me")):
            variants.append(manual_variant(entry, raw, [r_value, star, r_value, star], "Manual review: all * = H/Me; R = H/OMe/i-Pr.", f"R={r_value}_star={star}"))
        return with_variant_ids(entry, variants)

    if rid in {"LIT_0309", "LIT_0310", "LIT_0311"}:
        return with_variant_ids(entry, same_value_variants(entry, ["H", "OMe"], "Manual review: * = H/OMe."))

    return None


def apply_manual_review_overrides(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    curated: list[dict[str, Any]] = []
    manual_counts: Counter[str] = Counter()
    failures: list[dict[str, str]] = []
    for entry in entries:
        try:
            variants = manual_curated_variants(entry)
        except Exception as exc:
            variants = None
            failures.append(
                {
                    "reaction_id": entry["reaction_id"],
                    "raw_reaction_smiles": entry["raw_reaction_smiles"],
                    "error": str(exc),
                }
            )
        if variants is None:
            curated.append(entry)
            continue
        manual_counts[entry["reaction_id"]] = len(variants)
        curated.extend(variants)

    if failures:
        failure_path = OUT_DIR / "manual_review_override_failures.csv"
        with failure_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["reaction_id", "error", "raw_reaction_smiles"])
            writer.writeheader()
            writer.writerows(failures)

    summary = {
        "manual_override_parent_reactions": len(manual_counts),
        "manual_override_generated_entries": sum(manual_counts.values()),
        "manual_override_failures": failures,
        "manual_override_counts": dict(manual_counts),
        "broad_class_left_for_review": [],
    }
    (OUT_DIR / "manual_review_override_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return curated


def build_entries() -> list[dict[str, Any]]:
    md_files = sorted(SOURCE_DIR.glob("*.md"))
    entries: list[dict[str, Any]] = []
    seen_by_source = set()
    for path in md_files:
        for raw in extract_reactions_from_file(path):
            source_key = (str(path), raw.start_line, raw.raw)
            if source_key in seen_by_source:
                continue
            seen_by_source.add(source_key)
            if ">>" not in raw.raw:
                continue
            converted = reaction_to_smarts(raw.raw)
            conditions = parse_conditions(raw)
            classification = classify_reaction(raw, converted, conditions)
            reactant_components = split_components(converted["reactants_smiles"])
            product_components = split_components(converted["products_smiles"])
            has_generic = any(
                c["has_wildcard_or_generic_group"] for c in reactant_components + product_components
            )
            generic_context = " ".join(raw.context_before + raw.context_after + conditions["condition_lines"])
            generic_groups = find_symbol_meanings(raw.raw, generic_context) if has_generic else []
            resolution_status = generic_resolution_status(generic_groups, has_generic)
            rid = f"LIT_{len(entries) + 1:04d}"
            entry = {
                "reaction_id": rid,
                "name": conditions["nearby_caption"] or " / ".join(raw.heading_path[-2:]) or classification,
                "source": {
                    "file": str(raw.source_file),
                    "file_name": raw.source_file.name,
                    "line": raw.start_line,
                    "extraction_method": raw.extraction_method,
                    "heading_path": raw.heading_path,
                },
                "raw_reaction_smiles": raw.raw,
                "reactants": reactant_components,
                "products": product_components,
                "reaction_smarts": {
                    "forward": converted["forward_smarts"],
                    "retro": converted["retro_template"],
                },
                "classification": classification,
                "conditions": conditions,
                "generic_groups": generic_groups,
                "validation": {
                    "component_conversion_ok": converted["component_conversion_ok"],
                    "component_errors": converted["component_errors"],
                    "forward": converted["forward_validation"],
                    "retro": converted["retro_validation"],
                },
                "flags": {
                    "has_generic_group": has_generic,
                    "generic_resolution_status": resolution_status,
                    "aizynthfinder_ready": not has_generic,
                    "has_multiple_products": len(product_components) > 1,
                    "has_multiple_reactants": len(reactant_components) > 1,
                    "has_structured_conditions": bool(conditions["condition_text"]),
                },
                "curation": {
                    "confidence": 0.0,
                    "needs_review": not converted["retro_validation"]["ok"] or has_generic,
                    "review_reasons": [],
                },
            }
            reasons = entry["curation"]["review_reasons"]
            if has_generic:
                reasons.append(f"contains wildcard/generic group ({resolution_status})")
            if not converted["retro_validation"]["ok"]:
                reasons.append(converted["retro_validation"]["message"])
            if not conditions["condition_text"]:
                reasons.append("no condition text found nearby")
            entry["curation"]["confidence"] = confidence_for(entry)
            entries.append(entry)
    return entries


def write_json(entries: list[dict[str, Any]]) -> Path:
    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_dir": str(SOURCE_DIR),
        "description": (
            "Curated reaction entries extracted from retro_synthesis/files. "
            "reaction_smarts.forward is the normalized literature direction; "
            "reaction_smarts.retro is the reversed template used by AiZynthFinder."
        ),
        "entries": entries,
    }
    path = OUT_DIR / "flavonoid_literature_reactions.json"
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return path


def template_rows(entries: list[dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    full_rows = []
    for entry in entries:
        if not entry["validation"]["retro"]["ok"]:
            continue
        retro = entry["reaction_smarts"]["retro"]
        full_rows.append(
            {
                "template_code": len(full_rows),
                "retro_template": retro,
                "template_hash": hashlib.sha256(retro.encode()).hexdigest(),
                "classification": entry["classification"],
                "library_occurence": 1,
                "reaction_id": entry["reaction_id"],
                "name": entry["name"],
                "source_file": entry["source"]["file_name"],
                "source_line": entry["source"]["line"],
                "confidence": entry["curation"]["confidence"],
                "needs_review": entry["curation"]["needs_review"],
                "has_generic_group": entry["flags"]["has_generic_group"],
                "generic_resolution_status": entry["flags"]["generic_resolution_status"],
                "aizynthfinder_ready": entry["flags"]["aizynthfinder_ready"],
                "conditions": entry["conditions"]["condition_text"],
                "raw_reaction_smiles": entry["raw_reaction_smiles"],
                "forward_smarts": entry["reaction_smarts"]["forward"],
            }
        )
    full_df = pd.DataFrame(full_rows)
    if full_df.empty:
        return full_df, full_df
    grouped_rows = []
    for retro, group in full_df.groupby("retro_template", sort=False):
        first = group.iloc[0].to_dict()
        grouped_rows.append(
            {
                "retro_template": retro,
                "template_hash": hashlib.sha256(retro.encode()).hexdigest(),
                "classification": first["classification"],
                "library_occurence": int(len(group)),
            }
        )
    aizynth_df = pd.DataFrame(grouped_rows)
    return full_df, aizynth_df


def write_template_exports(entries: list[dict[str, Any]]) -> dict[str, Path]:
    full_df, aizynth_df = template_rows(entries)
    paths: dict[str, Path] = {}
    metadata_path = OUT_DIR / "flavonoid_literature_templates_metadata.csv"
    hdf5_path = OUT_DIR / "flavonoid_literature_templates.hdf5"
    csv_path = OUT_DIR / "flavonoid_literature_templates.csv.gz"
    if not full_df.empty:
        full_df.to_csv(metadata_path, index=False)
        full_df.to_csv(csv_path, index=False, compression="gzip")
        aizynth_df.to_hdf(hdf5_path, key="table", mode="w", format="table", data_columns=True)
        paths.update({"metadata_csv": metadata_path, "csv_gz": csv_path, "hdf5": hdf5_path})

        ready_full_df = full_df[full_df["aizynthfinder_ready"]].copy().reset_index(drop=True)
        ready_aizynth_rows = []
        if not ready_full_df.empty:
            for retro, group in ready_full_df.groupby("retro_template", sort=False):
                first = group.iloc[0].to_dict()
                ready_aizynth_rows.append(
                    {
                        "retro_template": retro,
                        "template_hash": hashlib.sha256(retro.encode()).hexdigest(),
                        "classification": first["classification"],
                        "library_occurence": int(len(group)),
                    }
                )
        ready_aizynth_df = pd.DataFrame(ready_aizynth_rows)
        ready_hdf5 = OUT_DIR / "flavonoid_literature_templates_aizynthfinder_ready.hdf5"
        ready_csv = OUT_DIR / "flavonoid_literature_templates_aizynthfinder_ready.csv.gz"
        if not ready_aizynth_df.empty:
            ready_aizynth_df.to_hdf(ready_hdf5, key="table", mode="w", format="table", data_columns=True)
            ready_full_df.to_csv(ready_csv, index=False, compression="gzip")
            paths.update({"ready_hdf5": ready_hdf5, "ready_csv_gz": ready_csv})

        root_hdf5 = TEMPLATE_DIR / "flavonoid_templates.hdf5"
        root_csv = TEMPLATE_DIR / "flavonoid_templates.csv.gz"
        root_export_df = ready_aizynth_df if not ready_aizynth_df.empty else aizynth_df
        root_metadata_df = ready_full_df if not ready_full_df.empty else full_df
        root_export_df.to_hdf(root_hdf5, key="table", mode="w", format="table", data_columns=True)
        root_metadata_df.to_csv(root_csv, index=False, compression="gzip")
        paths.update({"root_hdf5": root_hdf5, "root_csv_gz": root_csv})

        bio_df = root_export_df[root_export_df["classification"].eq("biosynthesis")].reset_index(drop=True)
        if not bio_df.empty:
            BIO_DIR.mkdir(parents=True, exist_ok=True)
            bio_path = BIO_DIR / "flavonoid_biosynthesis_templates.hdf5"
            bio_df.to_hdf(bio_path, key="table", mode="w", format="table", data_columns=True)
            paths["bio_hdf5"] = bio_path

        structural_df = root_export_df[~root_export_df["classification"].eq("biosynthesis")].reset_index(drop=True)
        if not structural_df.empty:
            STRUCTURAL_DIR.mkdir(parents=True, exist_ok=True)
            structural_path = STRUCTURAL_DIR / "flavonoid_structural_templates.hdf5"
            structural_df.to_hdf(structural_path, key="table", mode="w", format="table", data_columns=True)
            paths["structural_hdf5"] = structural_path
    return paths


def write_invalid(entries: list[dict[str, Any]]) -> Path:
    path = OUT_DIR / "flavonoid_literature_invalid_reactions.csv"
    rows = []
    for entry in entries:
        if entry["validation"]["retro"]["ok"]:
            continue
        rows.append(
            {
                "reaction_id": entry["reaction_id"],
                "source_file": entry["source"]["file_name"],
                "source_line": entry["source"]["line"],
                "classification": entry["classification"],
                "error": entry["validation"]["retro"]["message"],
                "component_errors": " | ".join(entry["validation"]["component_errors"]),
                "raw_reaction_smiles": entry["raw_reaction_smiles"],
            }
        )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "reaction_id",
                "source_file",
                "source_line",
                "classification",
                "error",
                "component_errors",
                "raw_reaction_smiles",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_generic_review(entries: list[dict[str, Any]]) -> dict[str, Path]:
    generic_entries = [entry for entry in entries if entry["flags"]["has_generic_group"]]
    json_path = OUT_DIR / "generic_group_review.json"
    csv_path = OUT_DIR / "generic_group_review.csv"
    md_path = OUT_DIR / "generic_group_review.md"

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "description": (
            "Manual-context review of placeholders from the markdown corpus. "
            "Entries are not placed in the default AiZynthFinder HDF5 while "
            "their reaction SMARTS still contains bare wildcard atoms."
        ),
        "entries": [
            {
                "reaction_id": entry["reaction_id"],
                "source": entry["source"],
                "classification": entry["classification"],
                "raw_reaction_smiles": entry["raw_reaction_smiles"],
                "generic_resolution_status": entry["flags"]["generic_resolution_status"],
                "aizynthfinder_ready": entry["flags"]["aizynthfinder_ready"],
                "generic_groups": entry["generic_groups"],
                "context_before": entry["conditions"]["context_before"],
                "context_after": entry["conditions"]["context_after"],
                "condition_text": entry["conditions"]["condition_text"],
            }
            for entry in generic_entries
        ],
    }
    json_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    fieldnames = [
        "reaction_id",
        "source_file",
        "line",
        "classification",
        "generic_resolution_status",
        "symbols",
        "meanings",
        "evidence",
        "raw_reaction_smiles",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for entry in generic_entries:
            annotations = entry["generic_groups"]
            writer.writerow(
                {
                    "reaction_id": entry["reaction_id"],
                    "source_file": entry["source"]["file_name"],
                    "line": entry["source"]["line"],
                    "classification": entry["classification"],
                    "generic_resolution_status": entry["flags"]["generic_resolution_status"],
                    "symbols": ";".join(ann["symbol"] for ann in annotations),
                    "meanings": ";".join(",".join(ann["meanings"]) for ann in annotations),
                    "evidence": " | ".join(ann["evidence"] for ann in annotations),
                    "raw_reaction_smiles": entry["raw_reaction_smiles"],
                }
            )

    status_counts = Counter(entry["flags"]["generic_resolution_status"] for entry in generic_entries)
    lines = [
        "# Generic Group Review",
        "",
        "This file records the context review for `[R]`, `[Ar]`, `[X]`, `[Y]`, and `*` placeholders.",
        "Entries in this list are retained in the full curated JSON/HDF5, but excluded from the default AiZynthFinder-ready export until the wildcard SMARTS is explicitly expanded or constrained.",
        "",
        "## Status Counts",
        "",
    ]
    for status, count in status_counts.most_common():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Reviewed Entries", ""])
    for entry in generic_entries:
        lines.append(
            f"### {entry['reaction_id']} - {entry['source']['file_name']}:{entry['source']['line']} "
            f"({entry['flags']['generic_resolution_status']})"
        )
        lines.append("")
        lines.append(f"- Classification: `{entry['classification']}`")
        lines.append(f"- Heading: {' / '.join(entry['source']['heading_path'])}")
        lines.append(f"- Raw: `{entry['raw_reaction_smiles']}`")
        if entry["generic_groups"]:
            for ann in entry["generic_groups"]:
                meanings = ", ".join(ann["meanings"])
                structures = ", ".join(
                    f"{st.get('label')}[{st.get('type')}]" for st in ann["possible_structures"]
                )
                lines.append(f"- `{ann['symbol']}` {ann.get('site') or ''}: {meanings}")
                lines.append(f"  Evidence: {ann['evidence']}")
                lines.append(f"  Structure: {structures}")
        else:
            lines.append("- No explicit nearby meaning could be extracted.")
        context = " ".join(entry["conditions"]["context_before"] + entry["conditions"]["context_after"])
        lines.append(f"- Context: {context[:800]}")
        lines.append("")
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"generic_review_json": json_path, "generic_review_csv": csv_path, "generic_review_md": md_path}


def write_stock(entries: list[dict[str, Any]]) -> dict[str, Path]:
    product_keys = {
        comp["inchikey"]
        for entry in entries
        for comp in entry["products"]
        if comp["inchikey"]
    }
    observed: dict[str, dict[str, Any]] = {}
    for entry in entries:
        for role, components in (("reactant", entry["reactants"]), ("product", entry["products"])):
            for comp in components:
                key = comp["inchikey"]
                if not key or comp["has_wildcard_or_generic_group"]:
                    continue
                row = observed.setdefault(
                    key,
                    {
                        "inchikey": key,
                        "smiles": comp["canonical_smiles"],
                        "mw": comp["molecular_weight"],
                        "observed_as_reactant": False,
                        "observed_as_product": False,
                        "source_reaction_ids": set(),
                        "source_files": set(),
                        "classifications": set(),
                    },
                )
                row[f"observed_as_{role}"] = True
                row["source_reaction_ids"].add(entry["reaction_id"])
                row["source_files"].add(entry["source"]["file_name"])
                row["classifications"].add(entry["classification"])

    rows = []
    strict_keys = []
    all_reactant_keys = []
    for key, row in sorted(observed.items(), key=lambda kv: kv[1]["smiles"]):
        is_initial_candidate = row["observed_as_reactant"] and key not in product_keys
        if row["observed_as_reactant"]:
            all_reactant_keys.append(key)
        if is_initial_candidate:
            strict_keys.append(key)
        rows.append(
            {
                "inchikey": key,
                "smiles": row["smiles"],
                "mw": row["mw"],
                "observed_as_reactant": row["observed_as_reactant"],
                "observed_as_product": row["observed_as_product"],
                "strict_starting_material_candidate": is_initial_candidate,
                "source_reaction_ids": ";".join(sorted(row["source_reaction_ids"])),
                "source_files": ";".join(sorted(row["source_files"])),
                "classifications": ";".join(sorted(row["classifications"])),
            }
        )

    candidates_path = OUT_DIR / "flavonoid_literature_stock_candidates.csv"
    fieldnames = [
        "inchikey",
        "smiles",
        "mw",
        "observed_as_reactant",
        "observed_as_product",
        "strict_starting_material_candidate",
        "source_reaction_ids",
        "source_files",
        "classifications",
    ]
    with candidates_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    strict_path = OUT_DIR / "flavonoid_literature_stock_inchikeys.txt"
    strict_path.write_text("\n".join(sorted(set(strict_keys))) + "\n", encoding="utf-8")
    all_reactants_path = OUT_DIR / "flavonoid_literature_all_reactants_stock_inchikeys.txt"
    all_reactants_path.write_text("\n".join(sorted(set(all_reactant_keys))) + "\n", encoding="utf-8")

    root_stock = TEMPLATE_DIR / "flavonoid_stock_inchikeys.txt"
    root_stock.write_text(strict_path.read_text(encoding="utf-8"), encoding="utf-8")

    if BIO_DIR.exists():
        bio_keys = {
            row["inchikey"]
            for row in rows
            if "biosynthesis" in row["classifications"] and row["strict_starting_material_candidate"]
        }
        (BIO_DIR / "flavonoid_biosynthesis_stock_inchikeys.txt").write_text(
            "\n".join(sorted(bio_keys)) + "\n", encoding="utf-8"
        )
    if STRUCTURAL_DIR.exists():
        structural_keys = {
            row["inchikey"]
            for row in rows
            if "biosynthesis" not in row["classifications"] and row["strict_starting_material_candidate"]
        }
        (STRUCTURAL_DIR / "flavonoid_structural_stock_inchikeys.txt").write_text(
            "\n".join(sorted(structural_keys)) + "\n", encoding="utf-8"
        )

    return {
        "candidates_csv": candidates_path,
        "strict_inchikeys": strict_path,
        "all_reactants_inchikeys": all_reactants_path,
        "root_stock": root_stock,
    }


def write_report(entries: list[dict[str, Any]], paths: dict[str, Path]) -> Path:
    classification_counts = Counter(entry["classification"] for entry in entries)
    source_counts = Counter(entry["source"]["file_name"] for entry in entries)
    valid_count = sum(1 for entry in entries if entry["validation"]["retro"]["ok"])
    generic_count = sum(1 for entry in entries if entry["flags"]["has_generic_group"])
    ready_count = sum(1 for entry in entries if entry["flags"].get("aizynthfinder_ready"))
    condition_count = sum(1 for entry in entries if entry["conditions"]["condition_text"])
    generic_status_counts = Counter(entry["flags"]["generic_resolution_status"] for entry in entries if entry["flags"]["has_generic_group"])
    hdf5_count = 0
    hdf5_path = paths.get("hdf5")
    if hdf5_path and hdf5_path.exists():
        hdf5_count = len(pd.read_hdf(hdf5_path, "table"))
    ready_hdf5_count = 0
    ready_hdf5_path = paths.get("ready_hdf5")
    if ready_hdf5_path and ready_hdf5_path.exists():
        ready_hdf5_count = len(pd.read_hdf(ready_hdf5_path, "table"))

    lines = [
        "# Flavonoid Literature Template Build Report",
        "",
        f"- Generated at: {datetime.now(timezone.utc).isoformat()}",
        f"- Source directory: `{SOURCE_DIR}`",
        f"- Extracted reaction entries: {len(entries)}",
        f"- RDKit-valid retro templates: {valid_count}",
        f"- Unique full templates: {hdf5_count}",
        f"- Conservative AiZynthFinder-ready templates: {ready_hdf5_count}",
        f"- Entries marked AiZynthFinder-ready: {ready_count}",
        f"- Entries with nearby condition text: {condition_count}",
        f"- Entries containing generic/wildcard groups: {generic_count}",
        "",
        "## Classification Counts",
        "",
    ]
    for key, count in classification_counts.most_common():
        lines.append(f"- {key}: {count}")
    lines.extend(["", "## Generic Resolution Counts", ""])
    if generic_status_counts:
        for key, count in generic_status_counts.most_common():
            lines.append(f"- {key}: {count}")
    else:
        lines.append("- none: 0")
    lines.extend(["", "## Source Counts", ""])
    for key, count in source_counts.most_common():
        lines.append(f"- {key}: {count}")
    lines.extend(["", "## Outputs", ""])
    for label, path in sorted(paths.items()):
        lines.append(f"- {label}: `{path}`")
    lines.extend(
        [
            "",
            "## Curation Notes",
            "",
            "- Manual review overrides from `generic_group_review.md` are applied before export: finite placeholders are expanded to concrete structures, enzyme/cofactor placeholders are moved to conditions, and broad aryl/alkyl placeholders are collapsed to parent skeleton templates.",
            "- The default `/home/ljx/retro_synthesis/templates/flavonoid_templates.hdf5` now includes every RDKit-valid curated entry; no generic placeholders remain in the ready export.",
            "- Full reviewed entries, including placeholder annotations, manual review notes, and possible structures, remain in `literature_curated/flavonoid_literature_reactions.json`; `generic_group_review.*` is regenerated as an empty residual review set when no generic entries remain.",
            "- `literature_curated/manual_review_override_summary.json` records how many variants each manually reviewed parent reaction generated.",
            "- One truncated source SMILES (`[I]c1ccc(OC(=O)C)`) is repaired to `[I]c1ccc(OC(=O)C)cc1` during extraction so the corresponding chalcone step remains usable.",
            "- The strict stock file includes valid reactant molecules that do not also appear as products in the extracted reaction set, to avoid marking route intermediates as purchasable by default.",
            "- `flavonoid_literature_reactions.json` keeps source context and extracted condition text for manual audit.",
        ]
    )
    report_path = OUT_DIR / "extraction_audit_report.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TEMPLATE_DIR.mkdir(parents=True, exist_ok=True)
    entries = apply_manual_review_overrides(build_entries())
    paths: dict[str, Path] = {}
    paths["curated_json"] = write_json(entries)
    summary_path = OUT_DIR / "manual_review_override_summary.json"
    if summary_path.exists():
        paths["manual_review_override_summary"] = summary_path
    paths.update(write_template_exports(entries))
    paths["invalid_csv"] = write_invalid(entries)
    paths.update(write_generic_review(entries))
    paths.update(write_stock(entries))
    paths["report"] = write_report(entries, paths)

    print(f"Extracted entries: {len(entries)}")
    print(f"Valid retro templates: {sum(1 for e in entries if e['validation']['retro']['ok'])}")
    for label, path in sorted(paths.items()):
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
