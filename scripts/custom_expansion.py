"""
Custom AiZynthFinder expansion strategy that applies a template library
uniformly, without a neural network policy model.

Usage in config.yml:

expansion:
  flavonoid:
    type: custom_expansion.UniformTemplateExpansion
    template: /home/ljx/retro_synthesis/templates/flavonoid_templates.hdf5
    use_rdchiral: false
    cutoff_number: 50
    min_atom_ratio: 0.8  # optional: pre-filter templates with poor atom balance
    min_retained_map_ratio: 0.8  # optional: post-filter by atom-map retention
"""
import re
from typing import List, Optional, Sequence, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from aizynthfinder.chem.reaction import RetroReaction, TemplatedRetroReaction
from aizynthfinder.context.policy.expansion_strategies import ExpansionStrategy


def _is_active_expansion(value) -> bool:
    """Interpret active_expansion values loaded from CSV/HDF5."""
    if pd.isna(value):
        return True
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value != 0
    if isinstance(value, str):
        return value.strip().lower() not in {"false", "0", "no", "n", "off"}
    return bool(value)


def _count_atoms(smarts_side: str) -> int:
    """Count heavy atoms in a SMARTS/SMILES string side.

    Handles both bracket expressions ([#6], [C@@H], [c], etc.)
    and non-bracket atoms (C, c, N, n, O, o, S, Cl, Br, etc.).
    """
    count = 0
    i = 0
    s = smarts_side
    while i < len(s):
        if s[i] == '[':
            # Bracket atom expression — count as 1 atom
            end = s.find(']', i)
            if end != -1:
                count += 1
                i = end + 1
            else:
                i += 1
        elif s[i].isalpha() and s[i] not in 'Hh':
            if s[i].isupper():
                # Handle two-letter elements
                if i + 1 < len(s) and s[i:i+2] in ('Cl', 'Br', 'Si', 'Se', 'Te'):
                    count += 1
                    i += 2
                elif s[i] in 'BCNOSPF':
                    count += 1
                    i += 1
                else:
                    i += 1
            elif s[i] in 'cnos':
                # Aromatic atoms
                count += 1
                i += 1
            else:
                i += 1
        else:
            i += 1
    return count


def _compute_atom_ratio(retro_template: str) -> float:
    """Compute the atom count ratio for a retro template.

    Returns precursor_atoms / target_atoms.
    For a well-formed retro template, this should be close to 1.0.
    Templates with ratio < 0.5 are severe "collapse" templates that
    transform a large target into a tiny precursor and should be filtered.

    In a retro template "precursor>>target":
    - precursor side (left) = what the target breaks into
    - target side (right) = what we're breaking from
    """
    if ">>" not in retro_template:
        return 0.0

    parts = retro_template.split(">>")
    if len(parts) != 2:
        return 0.0

    precursor_side, target_side = parts
    precursor_atoms = _count_atoms(precursor_side)
    target_atoms = _count_atoms(target_side)

    if target_atoms == 0:
        return 0.0

    return precursor_atoms / target_atoms


def _compute_map_retention(mapped_rxn_smiles: str) -> float:
    """Compute atom-map retention ratio for a mapped reaction SMILES.

    AiZynthFinder stores retro reactions as:
      current target molecule >> retrosynthetic precursors

    This computes:
      |target_map_numbers intersect precursor_map_numbers| / |target_map_numbers|

    A ratio of 0.0 means the reaction discards all mapped atoms from the target
    (scaffold collapse). A ratio >= 0.8 indicates good conservation.

    Returns 1.0 if no mapped atoms found (cannot audit, assume OK).
    """
    if ">>" not in mapped_rxn_smiles:
        return 1.0

    tgt_side, prec_side = mapped_rxn_smiles.split(">>", 1)
    tgt_maps = set(re.findall(r':(\d+)\]', tgt_side))
    prec_maps = set(re.findall(r':(\d+)\]', prec_side))

    if not tgt_maps:
        return 1.0  # No mapped atoms to audit

    return len(tgt_maps & prec_maps) / len(tgt_maps)


def _copy_without_atom_maps(mol: Chem.Mol) -> Chem.Mol:
    """Return a copy of an RDKit molecule with atom-map numbers cleared."""
    mol_copy = Chem.Mol(mol)
    for atom in mol_copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol_copy


class UniformTemplateExpansion(ExpansionStrategy):
    """
    Expansion strategy that returns every template whose reactant side matches
    the target molecule, all with the same prior probability.

    This is useful for evaluating a small, domain-specific template library
    when no trained policy model is available.
    """

    _required_kwargs = ["template"]

    def __init__(self, key: str, config, **kwargs):
        super().__init__(key, config, **kwargs)

        templatefile = kwargs["template"]
        if templatefile.endswith(".csv.gz") or templatefile.endswith(".csv"):
            self.templates: pd.DataFrame = pd.read_csv(
                templatefile, index_col=0, sep="\t"
            )
        else:
            self.templates = pd.read_hdf(templatefile, "table")

        if "active_expansion" in self.templates.columns:
            n_before_active = len(self.templates)
            active_mask = [
                _is_active_expansion(value)
                for value in self.templates["active_expansion"]
            ]
            self.templates = self.templates[active_mask].reset_index(drop=True)
            n_filtered_active = n_before_active - len(self.templates)
            if n_filtered_active > 0:
                self._logger.info(
                    f"Filtered {n_filtered_active}/{n_before_active} inactive "
                    "templates via active_expansion"
                )

        self.use_rdchiral: bool = bool(kwargs.get("use_rdchiral", False))
        self.cutoff_number: int = int(kwargs.get("cutoff_number", 50))
        # min_atom_ratio: pre-filter templates where precursor/target atom ratio is too low
        self.min_atom_ratio: float = float(
            kwargs.get("min_atom_ratio", kwargs.get("min_mapping_ratio", 0.0))
        )
        # min_retained_map_ratio: post-filter actions by atom-map retention
        # After applying the template, checks that target map numbers are preserved
        # in the precursor products. Rejects scaffold-collapse templates.
        self.min_retained_map_ratio: float = float(
            kwargs.get("min_retained_map_ratio", 0.0)
        )
        class_filter = kwargs.get("classifications", "")
        self._classifications: set = (
            {c.strip() for c in class_filter.split(",") if c.strip()}
            if isinstance(class_filter, str)
            else set(class_filter or [])
        )

        if self._classifications:
            self.templates = self.templates[
                self.templates["classification"].isin(self._classifications)
            ].reset_index(drop=True)

        # Pre-compute atom ratios and filter
        n_before = len(self.templates)
        if self.min_atom_ratio > 0:
            self._atom_ratios = [
                _compute_atom_ratio(row["retro_template"])
                for _, row in self.templates.iterrows()
            ]
            mask = [r >= self.min_atom_ratio for r in self._atom_ratios]
            self.templates = self.templates[mask].reset_index(drop=True)
            n_filtered = n_before - len(self.templates)
            if n_filtered > 0:
                self._logger.info(
                    f"Filtered {n_filtered}/{n_before} templates with "
                    f"atom ratio < {self.min_atom_ratio}"
                )
        else:
            self._atom_ratios = [1.0] * len(self.templates)

        self._logger.info(
            f"Loaded {len(self.templates)} custom templates for uniform expansion "
            f"(rdchiral={self.use_rdchiral}, cutoff={self.cutoff_number}, "
            f"min_atom_ratio={self.min_atom_ratio}, "
            f"min_retained_map_ratio={self.min_retained_map_ratio}, "
            f"classes={self._classifications or 'all'})"
        )

        self._rxns = [
            AllChem.ReactionFromSmarts(smarts)
            for smarts in self.templates["retro_template"]
        ]
        self._valid = [r is not None for r in self._rxns]

    def _check_map_retention(self, action: TemplatedRetroReaction) -> bool:
        """Check if an action preserves atom-map numbers above threshold.

        Applies the reaction and checks that the product's mapped atoms
        retain the target's map numbers. Returns True if retention >= threshold.

        If min_retained_map_ratio is 0, always returns True (disabled).
        """
        if self.min_retained_map_ratio <= 0:
            return True

        try:
            # Accessing reactants triggers AiZynthFinder's lazy reaction application.
            _ = action.reactants
            mapped_smiles = action.mapped_reaction_smiles()
            retention = _compute_map_retention(mapped_smiles)
            return retention >= self.min_retained_map_ratio
        except Exception:
            # If apply fails, reject the action
            return False

    def get_actions(
        self,
        molecules: Sequence,
        cache_molecules: Optional[Sequence] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        possible_actions: List[RetroReaction] = []
        priors: List[float] = []
        cache_molecules = cache_molecules or []

        for mol in molecules:
            n_added = 0
            n_rejected_map = 0
            unmapped_target = _copy_without_atom_maps(mol.mapped_mol)
            for idx, (move_index, move) in enumerate(self.templates.iterrows()):
                if not self._valid[idx]:
                    continue
                rxn = self._rxns[idx]
                if rxn.GetNumReactantTemplates() == 0:
                    continue

                reactant_template = _copy_without_atom_maps(rxn.GetReactantTemplate(0))
                if not unmapped_target.HasSubstructMatch(reactant_template):
                    continue

                metadata = dict(move)
                metadata.pop("retro_template", None)
                metadata["policy_probability"] = 1.0
                metadata["policy_probability_rank"] = n_added
                metadata["policy_name"] = self.key
                metadata["template_code"] = move_index
                metadata["template"] = move["retro_template"]

                action = TemplatedRetroReaction(
                    mol,
                    smarts=move["retro_template"],
                    metadata=metadata,
                    use_rdchiral=self.use_rdchiral,
                )

                # Post-filter: check atom-map retention
                if not self._check_map_retention(action):
                    n_rejected_map += 1
                    continue

                possible_actions.append(action)
                priors.append(1.0)
                n_added += 1
                if n_added >= self.cutoff_number:
                    break

            if n_rejected_map > 0:
                self._logger.info(
                    f"Rejected {n_rejected_map} actions for mol due to "
                    f"map retention < {self.min_retained_map_ratio}"
                )

        return possible_actions, priors
