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


# TODO More examples in docs for all classes
class ChemicalReaction:
    """
    A class that represents a chemical reaction and do operations on it.

    There are three calculation modes:

    1. The "force" mode is used when a user enters coefficients
    in the reaction string and wants the masses to be calculated
    whether the reaction is balanced or not.

    2. "check" mode is the same as force, but with reaction
    balance checks.

    3. "balance" mode  tries to automatically calculate
    coefficients from the reaction string.

    Important:
        Unlike other properties of this class, the [coefficients][chemsynthcalc.chemical_reaction.ChemicalReaction.coefficients]
        property can be set directly.

    Arguments:
        reaction (str): A reaction string
        mode (str): Coefficients calculation mode
        target (int): Index of target compound (0 by default, or first compound in the products), can be negative (limited by reactant)
        target_mass (float): Desired mass of target compound (in grams)
        precision (int): Value of rounding precision (8 by default)
        intify (bool): Is it required to convert the coefficients to integer values?

    Attributes:
        algorithm (str): Currently used calculation algorithm

    Raises:
        ValueError if precision or target mass <= 0
    """

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
            self.initial_reaction = reaction.replace(" ", "")

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
        return f"ChemicalReaction({self.reaction}, {self.mode}, {self.initial_target}, {self.target_mass}, {self.precision}, {self.intify})"

    def __str__(self) -> str:
        return self.reaction

    @property
    @lru_cache(maxsize=1)
    def reaction(self) -> str:
        """
        A string of chemical reaction.
        It is made a property to be relatively immutable.

        Returns:
            The reaction string
        """
        return self.initial_reaction

    @property
    @lru_cache(maxsize=1)
    def decomposed_reaction(self) -> ReactionDecomposer:
        """
        Decomposition of chemical reaction string and extraction of
        reaction separator, reactants, products and initial coefficients.

        Returns:
            A ReactionDecomposer object
        """
        return ReactionDecomposer(self.reaction)

    @property
    @lru_cache(maxsize=1)
    def _calculated_target(self) -> int:
        """
        Checks if initial_target is in the reaction's compounds range,
        and calculates the usable target integer.

        Returns:
            Final target

        Raises:
            IndexError if The target integer is not in the range
        """
        high = len(self.decomposed_reaction.products) - 1
        low = -len(self.decomposed_reaction.reactants)
        if self.initial_target <= high and self.initial_target >= low:
            return self.initial_target - low
        else:
            raise IndexError(
                f"The target integer {self.initial_target} should be in range {low} : {high}"
            )

    @property
    @lru_cache(maxsize=1)
    def chemformula_objs(self) -> list[ChemicalFormula]:
        """Decomposition of a list of formulas from the decomposed_reaction.

        Returns:
            Every compound as ChemicalFormula object
        """
        return [
            ChemicalFormula(formula, self.precision)
            for formula in self.decomposed_reaction.compounds
        ]

    @property
    @lru_cache(maxsize=1)
    def parsed_formulas(self) -> list[dict[str, float]]:
        """
        List of formulas parsed by [ChemicalFormulaParser][chemsynthcalc.formula_parser.ChemicalFormulaParser]
        """
        return [compound.parsed_formula for compound in self.chemformula_objs]

    @property
    @lru_cache(maxsize=1)
    def matrix(self) -> npt.NDArray[np.float64]:
        """Chemical reaction matrix.

        The first implementation of reaction matrix method is probably
        belongs to [Blakley](https://doi.org/10.1021/ed059p728). In general,
        a chemical reaction matrix is composed of the coefficients of each
        atom in each compound, giving a 2D array. The matrix composes
        naturally from previously parsed formulas.

        Returns:
            2D array of each atom amount in each formula
        """
        return ChemicalReactionMatrix(self.parsed_formulas).matrix

    @property
    @lru_cache(maxsize=1)
    def balancer(self) -> Balancer:
        """
        A balancer to  automatically balance chemical reaction by different matrix methods.

        Returns:
            A Balancer object
        """
        return Balancer(
            self.matrix,
            len(self.decomposed_reaction.reactants),
            self.precision,
            intify=self.intify,
        )

    @property
    @lru_cache(maxsize=1)
    def molar_masses(self) -> list[float]:
        """
        List of molar masses (in g/mol)

        Returns:
            List of molar masses of each compound in [chemformula_objs][chemsynthcalc.chemical_reaction.ChemicalReaction.chemformula_objs]"
        """
        return [compound.molar_mass for compound in self.chemformula_objs]

    @cached_property
    def coefficients(self) -> list[float | int] | list[int]:
        """
        Coefficients of the chemical reaction. Can be calculated,
        striped off the initial reaction string or set directly.

        Returns:
            A list of coefficients
        """
        coefs, self.algorithm = Coefficients(
            self.mode,
            self.parsed_formulas,
            self.matrix,
            self.balancer,
            self.decomposed_reaction,
        ).get_coefficients()
        return coefs

    @property
    @lru_cache(maxsize=1)
    def normalized_coefficients(self) -> list[float | int] | list[int]:
        """
        List of coefficients normalized on target compound.

        Target coefficient = 1.0

        Returns:
            Normalized coefficients
        """
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
        """
        Is the reaction balanced with the current coefficients?

        Returns:
            True if the reaction is balanced
        """
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
    @lru_cache(maxsize=1)
    def final_reaction(self) -> str:
        """
        Final representation of the reaction with coefficients.

        Returns:
            A string of the final reaction
        """
        return self._generate_final_reaction(self.coefficients)

    @property
    @lru_cache(maxsize=1)
    def final_reaction_normalized(self) -> str:
        """
        Final representation of the reaction with normalized coefficients.

        Returns:
            A string of the normalized final reaction
        """
        return self._generate_final_reaction(self.normalized_coefficients)

    @property
    @lru_cache(maxsize=1)
    def masses(self) -> list[float]:
        """
        List of masses of compounds (in grams).

        List of masses of the of formulas in reaction
        calculated with coefficients obtained by any of the 3 methods.
        Calculates masses by calculating amount of substance nu (nu=mass/molar mass).
        Coefficients of reaction are normalized to the target. After nu of target compound is
        calculated, it broadcasted to other compounds (with respect to their coefficients).

        Returns:
            A list of masses of compounds
        """
        nu = self.target_mass / self.molar_masses[self._calculated_target]
        masses = [
            round(molar * nu * self.normalized_coefficients[i], self.precision)
            for i, molar in enumerate(self.molar_masses)
        ]
        return masses

    @property
    @lru_cache(maxsize=1)
    def output_results(self) -> dict[str, object]:
        """
        Collection of every output of calculated ChemicalReaction properties.

        Returns:
            All outputs collected in one dictionary.
        """
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
        """
        Print a final result of calculations in stdout.

        Arguments:
            print_precision (int): print precision (4 digits by default)
        """
        ChemicalOutput(
            self.output_results, print_precision, obj=self.__class__.__name__
        ).print_results()

    def to_txt(self, filename: str = "default", print_precision: int = 4) -> None:
        """
        Export the final result of the calculations in a txt file.

        Arguments:
            filename (str): filename string (should end with .txt)
            print_precision (int): print precision (4 digits by default)
        """
        ChemicalOutput(
            self.output_results, print_precision, obj=self.__class__.__name__
        ).write_to_txt(filename)

    def to_json(self, print_precision: int = 4) -> str:
        """
        Serialization of output into JSON object.

        Arguments:
            print_precision (int): print precision (4 digits by default)

        Returns:
            A JSON-type object
        """
        return ChemicalOutput(
            self.output_results, print_precision, obj=self.__class__.__name__
        ).dump_to_json()

    def to_json_file(self, filename: str = "default", print_precision: int = 4) -> None:
        """
        Export a final result of calculations in a JSON file.

        Arguments:
            filename (str): filename string (should end with .json)
            print_precision (int): print precision (4 digits by default)
        """
        ChemicalOutput(
            self.output_results, print_precision, obj=self.__class__.__name__
        ).write_to_json_file(filename)
