import numpy as np
import numpy.typing as npt


class ChemicalReactionMatrix:
    def __init__(self, parsed_formulas: list[dict[str, float]]) -> None:
        self._parsed_formulas = parsed_formulas
        self._merged_dict: dict[str, float] = {
            k: v for d in self._parsed_formulas for k, v in d.items()
        }
        self._elements: list[str] = list(self._merged_dict.keys())
        self.matrix: npt.NDArray[np.float64] = self.create_reaction_matrix()

    def create_reaction_matrix(self) -> npt.NDArray[np.float64]:
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
