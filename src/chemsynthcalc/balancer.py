from fractions import Fraction

import numpy as np
import numpy.typing as npt

from .balancing_algos import BalancingAlgorithms
from .chem_errors import BalancingError


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
        self.coef_limit: int = 100000

    def _intify_coefficients(
        self, coefficients: list[float], max_denom: int = 100
    ) -> list[int]:
        ratios = np.array(
            [
                Fraction(val).limit_denominator(max_denom).as_integer_ratio()
                for val in coefficients
            ]
        )
        factor = np.lcm.reduce(ratios[:, 1])
        result = [round(v * factor) for v in coefficients]
        return result

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
            if np.allclose(reactants.sum(axis=0), products.sum(axis=0), rtol=tolerance):
                return True
            else:
                return False
        except Exception:
            return False

    def _calculate_by_method(self, method: str) -> list[float | int] | list[int]:
        match method:

            case "inv":
                coefficients: list[float] = np.round(
                    self._inv_algorithm(), decimals=self.round_precision
                ).tolist()

            case "gpinv":
                coefficients: list[float] = self._gpinv_algorithm().tolist()

            case "ppinv":
                coefficients: list[float] = self._ppinv_algorithm().tolist()

            case "comb":
                res: npt.NDArray[np.int32] | None = self._comb_algorithm()
                if res is not None:
                    return res.tolist()
                else:
                    raise BalancingError(f"Can't balance reaction by {method} method")

            case _:
                raise ValueError(f"No method {method}")

        if Balancer.is_reaction_balanced(
            self.reactant_matrix, self.product_matrix, coefficients
        ):
            if self.intify:
                intified = self._intify_coefficients(coefficients)
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

    '''
    def old_intify_coefficients(self, coefficients: list, limit: int) -> list:
        """Reduce the coefficients to integers by finding the greatest common divider.
        
        Arguments:
            coefficients (list): List of coefficients to intify
            limit (int): upper limit (max int coef)
        
        Returns:
            list: list of intified coefficients

        Example:
            >>> reaction = ChemicalReaction("P2O3+HClO3+H2O=H3PO4+HCl")
            >>> Balancer(reaction.reactant_matrix, reaction.product_matrix, 8, True, True).intify_coefficients([1.5, 1.0, 4.5, 3.0, 1.0], 100000)
            [3, 2, 9, 6, 2]
        """
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
        '''
