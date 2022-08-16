from .formula_parser import ChemicalFormulaParser
from .molar_mass import MolarMassCalculation

class ChemicalFormula():
    '''
    A base class for representing a single chemical formula. 
    It constructs with a formula string and calculates properties 
    parsed_formula, molar_mass, mass_percent, atomic_percent,
    oxide_percent from this string using ChemicalFormulaParser and
    MolarMassCalculation classes.
    '''

    def __init__(self, formula:str, rounding_order:int=5) -> None:
        self.formula:str = formula
        self.rounding_order:int = rounding_order

    def __repr__(self) -> str:
        return str(self.formula)

    def __str__(self) -> str:
        return str(self.formula)

    @property
    def parsed_formula(self) -> dict:
        '''
        Returns parsed dictionary representation of formula string
        created by ChemicalFormulaParser. For example, if formula is
        "K2SO4", then the parsed formula is {'S': 1.0, 'O': 4.0, 'K': 2.0}. 
        Order of atoms in dictionary is random due to the nature 
        of dictionaries in Python.
        '''
        return ChemicalFormulaParser(self.formula).parse_formula()

    @property
    def molar_mass(self) -> float:
        '''
        Molar mass of the formula, calculated from parsed formula
        using MolarMassCalculation. For example, molar mass of "K2SO4" 
        parsed into {'S': 1.0, 'O': 4.0, 'K': 2.0} is 174.2592 g/mol.
        '''
        return round(MolarMassCalculation(self.parsed_formula).calculate_molar_mass(), self.rounding_order)
    
    @property
    def mass_percent(self) -> dict:
        '''
        Calculates a mass percent (relative mass fraction)
        of atoms in parsed chemical formula. Return dictionary, for example
        output for "K2SO4" is {'K': 44.87373, 'O': 36.72552, 'S': 18.40075}.
        The values of mass content are in % (with 100% sum), not fraction
        '''
        output = MolarMassCalculation(self.parsed_formula).calculate_mass_percent()
        return {k: round(v, self.rounding_order) for k, v in output.items()}

    @property
    def atomic_percent(self) -> dict:
        '''
        Calculates an atomic percent (relative molar fraction) 
        dictionary of atoms in parsed chemical formula. For example, 
        output for "K2SO4" is {'S': 14.28571, 'K': 28.57143, 'O': 57.14286}.
        The values of mole content are in % (with 100% sum), not fraction
        '''
        output = MolarMassCalculation(self.parsed_formula).calculate_atomic_percent()
        return {k: round(v, self.rounding_order) for k, v in output.items()}

    @property
    def oxide_percent(self)  -> dict:
        '''
        Calculates an oxide percent (relative fraction of elements oxides) 
        dictionary of atoms in parsed chemical formula. Oxide types are listed
        in periodic_table file and can be change to any oxide formula. For example, 
        output for "K2SO4" is {'K2O': 54.05511, 'SO3': 45.94489}.
        The values of oxide content are in % (with 100% sum), not fraction
        '''
        output = MolarMassCalculation(self.parsed_formula).calculate_oxide_percent()
        return {k: round(v, self.rounding_order) for k, v in output.items()}