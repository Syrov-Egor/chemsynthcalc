import numpy as np
import numpy.typing as npt

from .balancer import Balancer
from .reaction_decomposer import ReactionDecomposer
from .chem_errors import (
    NoSuchMode,
    ReactionNotBalanced,
    BadCoeffiecients,
    ReactantProductDifference,
)
from .utils import to_integer


class Coefficients:
    """
    A class to calculate and validate reaction coefficients depending on the calculation "mode".

    Arguments:
        mode (str): Calculation mode ("force", "check", "balance")
        parsed_formulas (list[dict[str, float]]): List of formulas parsed by [ChemicalFormulaParser][chemsynthcalc.formula_parser.ChemicalFormulaParser]
        matrix (npt.NDArray[np.float64]): Reaction matrix created by [ChemicalReactionMatrix][chemsynthcalc.reaction_matrix.ChemicalReactionMatrix]
        balancer (Balancer): A [Balancer][chemsynthcalc.balancer.Balancer] object
        decomposed_reaction (ReactionDecomposer): A [ReactionDecomposer][chemsynthcalc.reaction_decomposer.ReactionDecomposer] object

    Attributes:
        initial_coefficients (list[float]): List of initial coefficients from decomposed reaction
    """

    def __init__(
        self,
        mode: str,
        parsed_formulas: list[dict[str, float]],
        matrix: npt.NDArray[np.float64],
        balancer: Balancer,
        decomposed_reaction: ReactionDecomposer,
    ) -> None:
        self.mode = mode
        self.matrix = matrix
        self.balancer = balancer
        self.initial_coefficients = decomposed_reaction.initial_coefficients
        self.decomposed_reaction = decomposed_reaction
        self.parsed_formulas = parsed_formulas

    def _calculate_coefficients(self) -> tuple[list[float | int] | list[int], str]:
        """
        Match a mode string and get coefficients depending on the mode

        Returns:
            Tuple of (coefficients, algorithm)

        Raise:
            [ReactionNotBalanced][chemsynthcalc.chem_errors.ReactionNotBalanced] if reaction is not balanced in the "check" mode <br />
            [NoSuchMode][chemsynthcalc.chem_errors.NoSuchMode] if there is no mode with that name
        """
        match self.mode:

            case "force":
                return (
                    to_integer(self.decomposed_reaction.initial_coefficients),
                    "user",
                )

            case "check":
                if Balancer.is_reaction_balanced(
                    self.balancer.reactant_matrix,
                    self.balancer.product_matrix,
                    self.decomposed_reaction.initial_coefficients,
                ):
                    return (
                        to_integer(self.decomposed_reaction.initial_coefficients),
                        "user",
                    )
                else:
                    raise ReactionNotBalanced("This reaction is not balanced!")

            case "balance":
                coefs, algorithm = self.balancer.auto()
                return (coefs, algorithm)

            case _:
                raise NoSuchMode(f"No mode {self.mode}")

    def coefficients_validation(
        self, coefficients: list[float | int] | list[int]
    ) -> None:
        """
        Validate a list of coefs.

        Arguments:
            coefficients (list[float | int] | list[int]): List of coefs

        Raise:
            [BadCoeffiecients][chemsynthcalc.chem_errors.BadCoeffiecients] if any coef <= 0 or
            lenght of list is not equal to the number of compounds.
        """
        if any(x <= 0 for x in coefficients):
            raise BadCoeffiecients("0 or -x in coefficients")
        elif len(coefficients) != self.matrix.shape[1]:
            raise BadCoeffiecients(
                f"Number of coefficients should be equal to {self.matrix.shape[1]}"
            )

    def _element_count_validation(self) -> None:
        """
        Calculate a [symmetric difference](https://en.wikipedia.org/wiki/Symmetric_difference)
        of two sets - left and right parts of the reaction. If this set is not empty, than
        some atoms are only in one part of the reaction (which is impossible).

        Raise:
            [ReactantProductDifference][chemsynthcalc.chem_errors.ReactantProductDifference] if diff set is not empty.
        """
        if self.mode != "force":
            reactants = {
                k: v
                for d in self.parsed_formulas[: len(self.decomposed_reaction.reactants)]
                for k, v in d.items()
            }
            products = {
                k: v
                for d in self.parsed_formulas[len(self.decomposed_reaction.reactants) :]
                for k, v in d.items()
            }
            diff = set(reactants.keys()) ^ (set(products.keys()))
            if diff:
                raise ReactantProductDifference(
                    f"Cannot balance this reaction, because element(s) {diff} are only in one part of the reaction"
                )

    def get_coefficients(self) -> tuple[list[float | int] | list[int], str]:
        """
        Validate atom's diff and coefs list and finally get a proper coefficients list.

        Returns:
            Tuple of (coefficients, algorithm)
        """
        self._element_count_validation()
        coefs = self._calculate_coefficients()
        self.coefficients_validation(coefs[0])
        return coefs
