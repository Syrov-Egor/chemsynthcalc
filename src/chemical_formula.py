from .formula_parser import ChemicalFormulaParser
from .molar_mass import MolarMassCalculation

class ChemicalFormula():
    """A base class for representing a single chemical formula. 
    It constructs with a formula string and calculates properties 
    parsed_formula, molar_mass and mass_percent from this string."""

    def __init__(self, formula:str, rounding_order:int=5) -> None:
        self.formula:str = formula
        self.rounding_order:int = rounding_order

    def __repr__(self):
        return str(self.formula)

    def __str__(self):
        return str(self.formula)

    @property
    def parsed_formula(self) -> dict:
        return ChemicalFormulaParser(self.formula).parse_formula()

    @property
    def molar_mass(self) -> float:
        return round(MolarMassCalculation(self.parsed_formula).calculate_molar_mass(), self.rounding_order)
    
    @property
    def mass_percent(self) -> dict:
        output = MolarMassCalculation(self.parsed_formula).calculate_mass_percent()
        return {k: round(v, self.rounding_order) for k, v in output.items()}

    @property
    def atomic_percent(self) -> dict:
        output = MolarMassCalculation(self.parsed_formula).calculate_atomic_percent()
        return {k: round(v, self.rounding_order) for k, v in output.items()}

    @property
    def oxide_percent(self)  -> dict:
        output = MolarMassCalculation(self.parsed_formula).calculate_oxide_percent()
        return {k: round(v, self.rounding_order) for k, v in output.items()}