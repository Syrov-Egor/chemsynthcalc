from functools import cached_property, lru_cache

import numpy as np
import numpy.typing as npt

from .reaction_validator import ReactionValidator
from .reaction_decomposer import ReactionDecomposer
from .chemical_formula import ChemicalFormula
from .reaction_matrix import ChemicalReactionMatrix
from .balancer import Balancer
from .coefs import Coefficients
from .chem_output import ChemicalOutput


# TODO Add a products-reagents diff function
class ChemicalReaction:
    def __init__(
        self,
        reaction: str = "",
        mode: str = "balance",
        target: int = 0,
        target_mass: float = 1.0,
        precision: int = 8,
        intify: bool = True,
    ) -> None:
        if ReactionValidator(reaction).validate_reaction():
            self.reaction = reaction.replace(" ", "")

        if precision > 0:
            self.precision: int = precision
        else:
            raise ValueError("precision <= 0")

        if target_mass > 0:
            self.target_mass: float = target_mass
        else:
            raise ValueError("target mass <= 0")

        self.intify: bool = intify
        self.mode: str = mode
        self.algorithm: str = "user"
        self.initial_target: int = target

    def __repr__(self) -> str:
        return "chemsynthcalc ChemicalReaction object: " + self.reaction

    def __str__(self) -> str:
        return self.reaction

    @property
    @lru_cache
    def decomposed_reaction(self) -> ReactionDecomposer:
        return ReactionDecomposer(self.reaction)

    @property
    @lru_cache
    def _calculated_target(self) -> int:
        high = len(self.decomposed_reaction.products) - 1
        low = -len(self.decomposed_reaction.reactants)
        if self.initial_target <= high and self.initial_target >= low:
            return self.initial_target + -low
        else:
            raise IndexError(
                f"The target integer {self.initial_target} should be in range {low} : {high}"
            )

    @property
    @lru_cache
    def chemformula_objs(self) -> list[ChemicalFormula]:
        return [
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

    @property
    @lru_cache
    def balancer(self) -> Balancer:
        return Balancer(
            self.matrix,
            len(self.decomposed_reaction.reactants),
            self.precision,
            intify=self.intify,
        )

    @property
    @lru_cache
    def molar_masses(self) -> list[float]:
        return [compound.molar_mass for compound in self.chemformula_objs]

    @cached_property
    def coefficients(self) -> list[float | int] | list[int]:
        coefs, self.algorithm = Coefficients(
            self.mode,
            self.parsed_formulas,
            self.matrix,
            self.balancer,
            self.decomposed_reaction,
        ).get_coefficients()
        return coefs

    @property
    @lru_cache
    def normalized_coefficients(self) -> list[float | int] | list[int]:
        target_compound = self.coefficients[self._calculated_target]
        normalized_coefficients: list[float | int] | list[int] = [
            coef / target_compound for coef in self.coefficients
        ]
        return [
            int(i) if i.is_integer() else round(i, self.precision)
            for i in normalized_coefficients
        ]

    @property
    def is_balanced(self) -> bool:
        return Balancer.is_reaction_balanced(
            self.balancer.reactant_matrix,
            self.balancer.product_matrix,
            self.coefficients,
        )

    def _generate_final_reaction(self, coefs: list[float | int] | list[int]) -> str:
        final_reaction = [
            (
                str(coefs[i]) + str(compound)
                if coefs[i] != 1 or coefs[i] != 1.0
                else str(compound)
            )
            for i, compound in enumerate(self.decomposed_reaction.compounds)
        ]
        final_reaction = (self.decomposed_reaction.reactant_separator).join(
            final_reaction
        )
        final_reaction = final_reaction.replace(
            self.decomposed_reaction.reactants[-1]
            + self.decomposed_reaction.reactant_separator,
            self.decomposed_reaction.reactants[-1] + self.decomposed_reaction.separator,
        )
        return final_reaction

    @property
    @lru_cache
    def final_reaction(self) -> str:
        return self._generate_final_reaction(self.coefficients)

    @property
    @lru_cache
    def final_reaction_normalized(self) -> str:
        return self._generate_final_reaction(self.normalized_coefficients)

    @property
    @lru_cache
    def masses(self) -> list[float]:
        nu = self.target_mass / self.molar_masses[self._calculated_target]
        masses = [
            round(molar * nu * self.normalized_coefficients[i], self.precision)
            for i, molar in enumerate(self.molar_masses)
        ]
        return masses

    @property
    @lru_cache
    def output_results(self) -> dict[str, object]:
        return {
            "initial reaction": self.reaction,
            "reaction matrix": self.matrix,
            "mode": self.mode,
            "formulas": self.decomposed_reaction.compounds,
            "coefficients": self.coefficients,
            "normalized coefficients": self.normalized_coefficients,
            "algorithm": self.algorithm,
            "is balanced": self.is_balanced,
            "final reaction": self.final_reaction,
            "final reaction normalized": self.final_reaction_normalized,
            "molar masses": self.molar_masses,
            "target": self.decomposed_reaction.compounds[self._calculated_target],
            "masses": self.masses,
        }

    def print_results(self, print_precision: int = 4) -> None:
        ChemicalOutput(
            self.output_results, print_precision, obj="reaction"
        ).print_results()

    def to_txt(self, filename: str = "default", print_precision: int = 4) -> None:
        ChemicalOutput(
            self.output_results, print_precision, obj="reaction"
        ).write_to_txt(filename)

    def to_json(self, print_precision: int = 4) -> str:
        return ChemicalOutput(
            self.output_results, print_precision, obj="reaction"
        ).dump_to_json()

    def to_json_file(self, filename: str = "default", print_precision: int = 4) -> None:
        ChemicalOutput(
            self.output_results, print_precision, obj="reaction"
        ).write_to_json_file(filename)
