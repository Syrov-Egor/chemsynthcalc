from functools import lru_cache

from .reaction_validator import ReactionValidator
from .reaction_decomposer import ReactionDecomposer
from .chemical_formula import ChemicalFormula

class ChemicalReaction:
    def __init__(self, reaction: str = "", precision: int = 8) -> None:
        if ReactionValidator(reaction).validate_reaction():
            self.reaction = reaction

        if precision > 0:
            self.precision: int = precision
        else:
            raise ValueError("precision <= 0")

        self.decomposed_reaction = ReactionDecomposer(self.reaction)
        self.chemformula_objs: list[ChemicalFormula] = [ChemicalFormula(formula, self.precision) for formula in self.decomposed_reaction.compounds]

    @property
    @lru_cache
    def parsed_formulas(self) -> list[dict[str, float]]:
        return [compound.parsed_formula for compound in self.chemformula_objs]
        