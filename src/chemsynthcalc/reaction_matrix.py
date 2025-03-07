import numpy as np
import numpy.typing as npt


class ChemicalReactionMatrix:
    """
    A class to create a dense float matrix from the parsed formulas.

    Arguments:
        parsed_formulas (list[dict[str, float]]): A list of formulas parsed by [ChemicalFormulaParser][chemsynthcalc.formula_parser.ChemicalFormulaParser]
    """

    def __init__(self, parsed_formulas: list[dict[str, float]]) -> None:
        self._parsed_formulas = parsed_formulas
        self._merged_dict: dict[str, float] = {
            k: v for d in self._parsed_formulas for k, v in d.items()
        }
        self._elements: list[str] = list(self._merged_dict.keys())
        self.matrix: npt.NDArray[np.float64] = self.create_reaction_matrix()

    def create_reaction_matrix(self) -> npt.NDArray[np.float64]:
        """
        Creates a 2D NumPy array from nested Python lists.
        The content of lists are exctracted from parsed dicts.

        Returns:
            A 2D NumPy array of the reaction matrix
        """
        matrix: list[list[float]] = []
        for element in self._elements:
            row: list[float] = []
            for compound in self._parsed_formulas:
                if element in compound.keys():
                    row.append(compound[element])
                else:
                    row.append(0.0)
            matrix.append(row)
        return np.array(matrix)
