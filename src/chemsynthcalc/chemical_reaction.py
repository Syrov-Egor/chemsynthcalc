from functools import lru_cache

import numpy as np
import numpy.typing as npt

from .reaction_validator import ReactionValidator
from .reaction_decomposer import ReactionDecomposer
from .chemical_formula import ChemicalFormula
from .reaction_matrix import ChemicalReactionMatrix


class ChemicalReaction:
    def __init__(self, reaction: str = "", precision: int = 8) -> None:
        if ReactionValidator(reaction).validate_reaction():
            self.reaction = reaction.replace(" ", "")

        if precision > 0:
            self.precision: int = precision
        else:
            raise ValueError("precision <= 0")

        self.decomposed_reaction = ReactionDecomposer(self.reaction)
        self.chemformula_objs: list[ChemicalFormula] = [
            ChemicalFormula(formula, self.precision)
            for formula in self.decomposed_reaction.compounds
        ]

    @property
    @lru_cache
    def parsed_formulas(self) -> list[dict[str, float]]:
        return [compound.parsed_formula for compound in self.chemformula_objs]

    @property
    @lru_cache
    def matrix(self) -> npt.NDArray[np.float64]:
        return ChemicalReactionMatrix(self.parsed_formulas).matrix
