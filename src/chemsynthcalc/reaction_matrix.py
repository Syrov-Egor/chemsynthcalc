from .periodic_table import periodic_table
import numpy as np


class ChemicalReactionMatrix:
    """
    Class that creates a reaction matrix from
    parsed formulas of the parsed reaction string.

    Parameters:
    * `parsed_all:list` - list of all parsed formulas

    For example:
    ```
    >>>parsed = ChemicalReaction("H2+O2=H2O").parsed_formulas
    >>>ChemicalReactionMatrix(parsed).create_reaction_matrix()
    [[2. 0. 2.]
     [0. 2. 1.]]
    ```
    """

    def __init__(self, parsed_all: list) -> None:
        self.parsed_all = parsed_all
        self.merged_dict = {k: v for d in self.parsed_all for k, v in d.items()}
        self.electronegativity = [(x[0], x[2]) for x in periodic_table]
        self.electronegativity_pattern = [
            atom[0] for atom in sorted(self.electronegativity, key=lambda tup: tup[1])
        ]
        self.list_of_atoms = [x[0] for x in periodic_table]

    def create_reaction_matrix(self) -> np.array:
        """
        Method that creates a chemical matrix of reaction from list of
        parsed compound. Initially contains a species type (atom label)
        at the start of each row of matrix (might be used in case of future
        expansion of matrix representation in print). Retruns 2D
        NumPy array (since Matrix class is no longer recommended
        https://numpy.org/doc/stable/reference/generated/numpy.matrix.html).
        """
        elements = list(self.merged_dict.keys())
        elements_vector = np.array(elements, dtype="str_").T
        matrix = []
        for i, element in enumerate(elements):
            row = []
            for parsed in self.parsed_all:
                if element in parsed.keys():
                    row.append(parsed.get(elements[i]))
                else:
                    row.append(0)
            matrix.append(row)
        return np.array(matrix)
