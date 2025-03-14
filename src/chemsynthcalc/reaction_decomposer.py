from .reaction import Reaction


class ReactionDecomposer(Reaction):
    """
    Decomposition of chemical reaction string and extraction of
    reaction separator, reactants, products and initial coefficients.

    Arguments:
        reaction (str): A reaction string

    Attributes:
        separator (str): A reactants - products separator (usually "+")
        initial_coefficients (list[float]): A list of coefficients striped from the formulas
        reactants (list[str]): A list of compound from letf part of the reaction equation
        products (list[str]): A list of compound from right part of the reaction equation
        compounds (list[str]): reactants + products
    """

    def __init__(self, reaction: str) -> None:
        super().__init__(reaction)

        self.separator: str = self.extract_separator()

        self._initial_reactants: list[str] = self.reaction.split(self.separator)[
            0
        ].split(self.reactant_separator)
        self._initial_products: list[str] = self.reaction.split(self.separator)[
            1
        ].split(self.reactant_separator)
        self._splitted_compounds: list[tuple[float, str]] = [
            self.split_coefficient_from_formula(formula)
            for formula in self._initial_reactants + self._initial_products
        ]

        self.initial_coefficients: list[float] = [
            atom[0] for atom in self._splitted_compounds
        ]
        self.compounds: list[str] = [atom[1] for atom in self._splitted_compounds]
        self.reactants: list[str] = self.compounds[: len(self._initial_reactants)]
        self.products: list[str] = self.compounds[len(self._initial_reactants) :]

    def __str__(self) -> str:
        return f"separator: {self.separator}; reactants: {self.reactants}; products: {self.products}"

    def __repr__(self) -> str:
        return f"ReactionDecomposer({self.reaction})"

    def split_coefficient_from_formula(self, formula: str) -> tuple[float, str]:
        """
        Split the coefficient (int or float) from string containing formula and coef.

        Parameters:
            formula (str): Formula string

        Returns:
            A tuple of (coefficient, formula)
        """
        if not formula[0].isdigit():
            return 1.0, formula
        else:
            coef: list[str] = []
            i: int = 0
            for i, symbol in enumerate(formula):
                if symbol.isdigit() or symbol == ".":
                    coef.append(symbol)
                else:
                    break
            return float("".join(coef)), formula[i:]
