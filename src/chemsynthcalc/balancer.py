from fractions import Fraction

import numpy as np
import numpy.typing as npt

from .balancing_algos import BalancingAlgorithms
from .chem_errors import BalancingError
from .utils import find_lcm, find_gcd


class Balancer(BalancingAlgorithms):
    def __init__(
        self,
        matrix: npt.NDArray[np.float64],
        separator_pos: int,
        round_precision: int,
        intify: bool = True,
    ) -> None:
        super().__init__(matrix, separator_pos)

        if round_precision > 0:
            self.round_precision: int = round_precision
        else:
            raise ValueError("precision <= 0")

        self.intify: bool = intify
        self.coef_limit: int = 1_000_000

    def _intify_coefficients(
        self, coefficients: list[float], limit: int
    ) -> list[float | int] | list[int]:
        initial_coefficients = coefficients
        frac = [Fraction(x).limit_denominator() for x in coefficients]
        vals = [
            int(
                fr.numerator
                * find_lcm([fr.denominator for fr in frac])
                / fr.denominator
            )
            for fr in frac
        ]
        coefficients = [int(val / find_gcd(vals)) for val in vals]
        if any(x > limit for x in coefficients):
            return initial_coefficients
        return coefficients

    @staticmethod
    def is_reaction_balanced(
        reactant_matrix: npt.NDArray[np.float64],
        product_matrix: npt.NDArray[np.float64],
        coefficients: list[float] | list[int],
        tolerance: float = 1e-8,
    ) -> bool:
        try:
            reactants = np.multiply(
                reactant_matrix.T,
                np.array(coefficients)[: reactant_matrix.shape[1], None],
            )
            products = np.multiply(
                product_matrix.T,
                np.array(coefficients)[reactant_matrix.shape[1] :, None],
            )
            return np.allclose(
                reactants.sum(axis=0), products.sum(axis=0), rtol=tolerance
            )

        except Exception:
            return False

    def _calculate_by_method(self, method: str) -> list[float | int] | list[int]:
        match method:

            case "inv":
                coefficients: list[float] = np.round(
                    self._inv_algorithm(), decimals=self.round_precision
                ).tolist()

            case "gpinv":
                coefficients: list[float] = np.round(
                    self._gpinv_algorithm(), decimals=self.round_precision + 2
                ).tolist()

            case "ppinv":
                coefficients: list[float] = np.round(
                    self._ppinv_algorithm(), decimals=self.round_precision + 2
                ).tolist()

            case "comb":
                res: npt.NDArray[np.int32] | None = self._comb_algorithm()
                if res is not None:
                    return res.tolist()
                else:
                    raise BalancingError(f"Can't balance reaction by {method} method")

            case _:
                raise ValueError(f"No method {method}")

        if (
            Balancer.is_reaction_balanced(
                self.reactant_matrix, self.product_matrix, coefficients
            )
            and all(x > 0 for x in coefficients)
            and len(coefficients) == self.reaction_matrix.shape[1]
        ):
            if self.intify:
                intified = self._intify_coefficients(coefficients, self.coef_limit)
                if all(x < self.coef_limit for x in intified):
                    return intified
                else:
                    return coefficients
            else:
                return coefficients
        else:
            raise BalancingError(f"Can't balance reaction by {method} method")

    def inv(self) -> list[float | int] | list[int]:
        return self._calculate_by_method("inv")

    def gpinv(self) -> list[float | int] | list[int]:
        return self._calculate_by_method("gpinv")

    def ppinv(self) -> list[float | int] | list[int]:
        return self._calculate_by_method("ppinv")

    def comb(self) -> list[float | int] | list[int]:
        return self._calculate_by_method("comb")

    def auto(self) -> tuple[list[float | int] | list[int], str]:
        try:
            return (self.inv(), "inverse")
        except Exception:
            pass
        try:
            return (self.gpinv(), "general pseudoinverse")
        except Exception:
            pass
        try:
            return (self.gpinv(), "partial pseudoinverse")
        except Exception:
            raise BalancingError("Can't balance this reaction by any method")
