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
"""
from typing import List, Optional, Sequence, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from aizynthfinder.chem.reaction import RetroReaction, TemplatedRetroReaction
from aizynthfinder.context.policy.expansion_strategies import ExpansionStrategy


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

        self.use_rdchiral: bool = bool(kwargs.get("use_rdchiral", False))
        self.cutoff_number: int = int(kwargs.get("cutoff_number", 50))
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

        self._logger.info(
            f"Loaded {len(self.templates)} custom templates for uniform expansion "
            f"(rdchiral={self.use_rdchiral}, cutoff={self.cutoff_number}, "
            f"classes={self._classifications or 'all'})"
        )

        self._rxns = [
            AllChem.ReactionFromSmarts(smarts)
            for smarts in self.templates["retro_template"]
        ]
        self._valid = [r is not None for r in self._rxns]

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
            for idx, (move_index, move) in enumerate(self.templates.iterrows()):
                if not self._valid[idx]:
                    continue
                rxn = self._rxns[idx]
                if rxn.GetNumReactantTemplates() == 0:
                    continue

                reactant_template = rxn.GetReactantTemplate(0)
                if not mol.mapped_mol.HasSubstructMatch(reactant_template):
                    continue

                metadata = dict(move)
                metadata.pop("retro_template", None)
                metadata["policy_probability"] = 1.0
                metadata["policy_probability_rank"] = n_added
                metadata["policy_name"] = self.key
                metadata["template_code"] = move_index
                metadata["template"] = move["retro_template"]

                possible_actions.append(
                    TemplatedRetroReaction(
                        mol,
                        smarts=move["retro_template"],
                        metadata=metadata,
                        use_rdchiral=self.use_rdchiral,
                    )
                )
                priors.append(1.0)
                n_added += 1
                if n_added >= self.cutoff_number:
                    break

        return possible_actions, priors
