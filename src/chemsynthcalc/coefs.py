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

    def _coefficients_validation(
        self, coefficients: list[float | int] | list[int]
    ) -> None:
        if any(x <= 0 for x in coefficients):
            raise BadCoeffiecients("0 or -x in coefficients")
        elif len(coefficients) != self.matrix.shape[1]:
            raise BadCoeffiecients(
                f"Number of coefficients should be equal to {self.matrix.shape[1]}"
            )

    def _element_count_validation(self) -> None:
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
        self._element_count_validation()
        coefs = self._calculate_coefficients()
        self._coefficients_validation(coefs[0])
        return coefs
