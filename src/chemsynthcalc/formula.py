class ChemicalFormula:

    def __init__(self, formula: str) -> None:
        self.atom_regex: str = r"([A-Z][a-z]*)"
        self.coefficient_regex: str = r"((\d+(\.\d+)?)*)"
        self.allowed_symbols: str = r"[^A-Za-z0-9.({[)}\]*·•]"
        self.atom_and_coefficient_regex: str = self.atom_regex + self.coefficient_regex
        self.opener_brackets: str = "({["
        self.closer_brackets: str = ")}]"
        self.adduct_symbols: str = "*·•"

        self.formula: str = formula
