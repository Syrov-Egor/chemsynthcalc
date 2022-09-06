from .formula_parser import ChemicalFormulaParser
from .molar_mass import MolarMassCalculation
from .chem_output import FormulaOutput
from .chemutils import arguments_type_checking
from functools import lru_cache

class ChemicalFormula():
    '''
    A base class for representing a single chemical formula. 
    It constructs with a formula string and calculates properties 
    parsed_formula, molar_mass, mass_percent, atomic_percent,
    oxide_percent from this string using ChemicalFormulaParser and
    MolarMassCalculation classes.

    Parameters:
    * `formula:srt` - string of chemical formula
    * `rounding_order:int` - value of rounding order (5 by default)
    '''

    def __init__(self, 
    formula:str = "", 
    rounding_order:int=8) -> None:
        arguments_type_checking(formula, str)
        arguments_type_checking(rounding_order, int)
        if formula == "":
            raise ValueError("No formula!")
        self.formula:str = formula.replace(" ", "")
        if rounding_order > 0:
            self.rounding_order:int = rounding_order
        else:
            raise ValueError("rounding order <= 0")

    def __repr__(self) -> str:
        return str(self.formula)

    def __str__(self) -> str:
        return str(self.formula)

    @property
    @lru_cache
    def parsed_formula(self) -> dict:
        '''
        Returns parsed dictionary representation of formula string
        created by ChemicalFormulaParser. For example, if formula is
        "K2SO4", then the parsed formula is {'K': 2.0, 'S': 1.0, 'O': 4.0}. 
        '''
        parsed = ChemicalFormulaParser(self.formula).parse_formula()
        return {k: round(v, self.rounding_order+3) for k, v in parsed.items()}
    
    @property
    @lru_cache
    def molar_mass(self) -> float:
        '''
        Molar mass of the formula, calculated from parsed formula
        using MolarMassCalculation. For example, molar mass of "K2SO4" 
        parsed into {'K': 2.0, 'S': 1.0, 'O': 4.0} is 174.2592 g/mol.
        '''
        return round(MolarMassCalculation(self.parsed_formula).calculate_molar_mass(), self.rounding_order)
    
    @property
    @lru_cache
    def mass_percent(self) -> dict:
        '''
        Calculates a mass percent (relative mass fraction)
        of atoms in parsed chemical formula. Return dictionary, for example
        output for "K2SO4" is {'K': 44.87373, 'S': 18.40075, 'O': 36.72552}.
        The values of mass content are in % (with 100% sum), not fraction
        '''
        output = MolarMassCalculation(self.parsed_formula).calculate_mass_percent()
        return {k: round(v, self.rounding_order) for k, v in output.items()}

    @property
    @lru_cache
    def atomic_percent(self) -> dict:
        '''
        Calculates an atomic percent (relative molar fraction) 
        dictionary of atoms in parsed chemical formula. For example, 
        output for "K2SO4" is {'K': 28.57143, 'S': 14.28571, 'O': 57.14286}.
        The values of mole content are in % (with 100% sum), not fraction
        '''
        output = MolarMassCalculation(self.parsed_formula).calculate_atomic_percent()
        return {k: round(v, self.rounding_order) for k, v in output.items()}

    @property
    @lru_cache
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
    
    @property
    @lru_cache
    def output_results(self) -> dict:
        '''
        Dictionary of calculation resulst output for
        `ChemicalFormula`.
        '''
        output:dict = {
            "formula" : self.formula,
            "parsed formula" : self.parsed_formula,
            "molar mass" : self.molar_mass,
            "mass percent" : self.mass_percent,
            "atomic percent" : self.atomic_percent,
            "oxide percent" : self.oxide_percent
        }
        return output

    def print_results(self, print_rounding_order:int = 4) -> None:
        '''
        Method to print a final result of calculations
        in terminal.
        '''
        printing = FormulaOutput(self.output_results).print_results(print_rounding_order)
        return
        
    def export_to_txt(self, filename:str='default', print_rounding_order:int = 4)  -> None:
        '''
        Method to print a final result of calculations
        in txt file.
        '''
        printing = FormulaOutput(self.output_results).export_to_txt(filename, print_rounding_order)
        return

    def as_json(self, print_rounding_order:int = 4) -> str:
        '''
        Serialization of output into JSON object
        '''
        return FormulaOutput(self.output_results).dump_to_json(print_rounding_order)
        
    def export_to_json(self, filename:str='default', print_rounding_order:int = 4)  -> None:
        '''
        Method to print a final result of calculations
        in JSON file.
        '''
        printing = FormulaOutput(self.output_results).export_to_json(filename, print_rounding_order)
        return