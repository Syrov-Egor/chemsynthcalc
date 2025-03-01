class Formula:
    """
    A base class for
    [ChemicalFormulaParser][chemsynthcalc.formula_parser.ChemicalFormulaParser] and
    [FormulaValidator][chemsynthcalc.formula_validator.FormulaValidator]
    containing regexes and symbols.

    Parameters:
        formula (str): Formula string

    Attributes:
        atom_regex (str): Regular expression for finding atoms in formula
        coefficient_regex (str): Regular expression for atoms amounts in formula
        atom_and_coefficient_regex (str): atom_regex+coefficient_regex
        opener_brackets (str): Opener brackets variations
        closer_brackets (str): Closer brackets variations
        adduct_symbols (str): Symbols for adduct notation (water of crystallization most often)
    """

    def __init__(self, formula: str) -> None:
        self.atom_regex: str = r"([A-Z][a-z]*)"
        self.coefficient_regex: str = r"((\d+(\.\d+)?)*)"
        self.allowed_symbols: str = r"[^A-Za-z0-9.({[)}\]*·•]"
        self.atom_and_coefficient_regex: str = self.atom_regex + self.coefficient_regex
        self.opener_brackets: str = "({["
        self.closer_brackets: str = ")}]"
        self.adduct_symbols: str = "*·•"

        self.formula: str = formula.replace(" ", "")
