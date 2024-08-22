from fractions import Fraction

import numpy as np
import numpy.typing as npt

from .balancing_algos import BalancingAlgorithms
#from .utils import find_gcd, find_lcm

class Balancer(BalancingAlgorithms):
    def __init__(self, matrix: npt.NDArray[np.float64], separator_pos: int, round_precision: int, intify: bool = True) -> None:
        super().__init__(matrix, separator_pos)

        if round_precision > 0:
            self.round_precision: int = round_precision
        else:
            raise ValueError("precision <= 0")
        
        self.intify: bool = intify

    def intify_coefficients(self, coefficients: list[float], max_denom: int = 1000) -> list[int]:
        ratios = np.array([Fraction(val).limit_denominator(max_denom).as_integer_ratio() for val in coefficients])
        factor = np.lcm.reduce(ratios[:,1])
        result = [round(v * factor) for v in coefficients]
        return result
    
    def inv(self) -> list[float] | list[int]:
        coefficients: list[float] = np.round(self._inv_algorithm(), decimals=self.round_precision).tolist()
        if self.intify:
            return self.intify_coefficients(coefficients)
        else:
            return coefficients
        
    def gpinv(self) -> list[float] | list[int]:
        coefficients: list[float] = self._gpinv_algorithm().tolist()
        if self.intify:
            return self.intify_coefficients(coefficients)
        else:
            return coefficients
        
    def ppinv(self) -> list[float] | list[int]:
        coefficients: list[float] = self._ppinv_algorithm().tolist()
        if self.intify:
            return self.intify_coefficients(coefficients)
        else:
            return coefficients

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