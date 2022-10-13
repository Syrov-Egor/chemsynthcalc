from .formula_parser import ChemicalFormulaParser
from .molar_mass import MolarMassCalculation
from .chem_output import FormulaOutput
from .chemutils import arguments_type_checking
from functools import lru_cache


class ChemicalFormula:
    """A base class for representing a single chemical formula.

    It constructs with a formula string and calculates properties:
    parsed formula, molar mass, mass percent, atomic percent,
    oxide percent  from this string using 
    :class:`chemsynthcalc.formula_parser.ChemicalFormulaParser` and
    :class:`chemsynthcalc.molar_mass.MolarMassCalculation`.

    Arguments:
        formula (srt): string of chemical formula
        rounding_order (int): value of rounding precision (8 by default)

    Raises:
        ValueError: if the input string is empty
        ValueError: if rounding order <= 0

    Examples:
        >>> ChemicalFormula("H2O")
        H2O
        >>> ChemicalFormula("H2O").molar_mass
        18.015
        >>> ChemicalFormula("H2O").mass_percent
        {'H': 11.19067444, 'O': 88.80932556}
    """

    def __init__(self, formula: str = "", rounding_order: int = 8) -> None:
        arguments_type_checking(formula, str)
        arguments_type_checking(rounding_order, int)
        if formula == "":
            raise ValueError("No formula!")
        self.formula: str = formula.replace(" ", "")
        if rounding_order > 0:
            self.rounding_order: int = rounding_order
        else:
            raise ValueError("rounding order <= 0")

    def __repr__(self) -> str:
        return str(self.formula)

    def __str__(self) -> str:
        return str(self.formula)

    @property
    @lru_cache
    def parsed_formula(self) -> dict:
        """Formula parsed into dict.

        Returns:
            dict: Parsed dictionary representation of formula string
            created by :class:`chemsynthcalc.formula_parser.ChemicalFormulaParser`.

        Example:
            >>> ChemicalFormula("K2SO4").parsed_formula
            {'K': 2.0, 'S': 1.0, 'O': 4.0}
        """
        parsed = ChemicalFormulaParser(self.formula).parse_formula()
        return {k: round(v, self.rounding_order + 3) for k, v in parsed.items()}

    @property
    @lru_cache
    def molar_mass(self) -> float:
        """Molar mass of formula.

        Returns:
            float: The `molar mass <https://en.wikipedia.org/wiki/Molar_mass>`_
            of the formula (in g/mol), calculated from
            parsed the formula using :class:`chemsynthcalc.molar_mass.MolarMassCalculation`.

        Example:
            >>> ChemicalFormula("K2SO4").molar_mass
            174.252
        """
        return round(
            MolarMassCalculation(self.parsed_formula).calculate_molar_mass(),
            self.rounding_order,
        )

    @property
    @lru_cache
    def mass_percent(self) -> dict:
        """The percentage of mass of atoms in the formula.

        Returns:
            dict: A mass percent or 
            `relative mass fraction <https://en.wikipedia.org/wiki/Mass_fraction_(chemistry)>`_
            of atoms in parsed chemical formula. The values of
            mass content are in % (with 100% sum), not fraction.
        
        Example:
            >>> ChemicalFormula("K2SO4").mass_percent
            {'K': 44.87523816, 'S': 18.39864105, 'O': 36.72612079}
        """
        output = MolarMassCalculation(self.parsed_formula).calculate_mass_percent()
        return {k: round(v, self.rounding_order) for k, v in output.items()}

    @property
    @lru_cache
    def atomic_percent(self) -> dict:
        """Atomic percents of atoms in the formula.

        Returns:
            dict: An atomic percent or 
            `relative mole fraction <https://en.wikipedia.org/wiki/Mole_fraction>`_
            dictionary of atoms in a parsed chemical formula.
            The values of mole content are in % (with 100% sum),
            not fraction.
        
        Example:
            >>> ChemicalFormula("K2SO4").atomic_percent
            {'K': 28.57142857, 'S': 14.28571429, 'O': 57.14285714}
        """
        output = MolarMassCalculation(self.parsed_formula).calculate_atomic_percent()
        return {k: round(v, self.rounding_order) for k, v in output.items()}

    @property
    @lru_cache
    def oxide_percent(self) -> dict:
        """Oxide percents of metals in formula.

        Returns:
            dict: An oxide percent or 
            `oxide fraction <https://d32ogoqmya1dw8.cloudfront.net/files/introgeo/studio/examples/minex02.pdf>`_ 
            dictionary of atoms in parsed chemical formula. Oxide types are listed
            in the :mod:`chemsynthcalc.periodic_table` file and can be changed to
            any oxide formula. The values of oxide content
            are in % (with 100% sum), not fraction.

        Example:
            >>> ChemicalFormula("K2SO4").oxide_percent
            {'K2O': 54.05676836, 'SO3': 45.94323164}
        """
        output = MolarMassCalculation(self.parsed_formula).calculate_oxide_percent()
        return {k: round(v, self.rounding_order) for k, v in output.items()}

    @property
    @lru_cache
    def output_results(self) -> dict:
        """Dictionary of the calculation result output for class.

        Returns:
            dict: Output dictionary for all properties listed above.

        Example:
            >>> ChemicalFormula("H2O").output_results
            {'formula': 'H2O', 
            'parsed formula': {'H': 2.0, 'O': 1.0}, 
            'molar mass': 18.015, 
            'mass percent': {'H': 11.19067444, 'O': 88.80932556}, 
            'atomic percent': {'H': 66.66666667, 'O': 33.33333333}, 
            'oxide percent': {'H2O': 100.0}}
        """
        output: dict = {
            "formula": self.formula,
            "parsed formula": self.parsed_formula,
            "molar mass": self.molar_mass,
            "mass percent": self.mass_percent,
            "atomic percent": self.atomic_percent,
            "oxide percent": self.oxide_percent,
        }
        return output

    def print_results(self, print_rounding_order: int = 4) -> None:
        """Method to print a final result of calculations in terminal.

        Arguments:
            print_rounding_order (int): print precision (4 digits by default)
                    
        Returns:
            None
        """
        printing = FormulaOutput(self.output_results).print_results(
            print_rounding_order
        )
        return

    def export_to_txt(
        self, filename: str = "default", print_rounding_order: int = 4
    ) -> None:
        """Method to print the final result of the calculations in a txt file.

        Arguments:
            filename (str): filename string (should end with .txt)
            print_rounding_order (int): print precision (4 digits by default)
        
        Returns:
            None
        """
        printing = FormulaOutput(self.output_results).export_to_txt(
            filename, print_rounding_order
        )
        return

    def as_json(self, print_rounding_order: int = 4) -> str:
        """Serialization of output into JSON object.

        Arguments:
            print_rounding_order (int): print precision (4 digits by default)
        
        Returns:
            str: A JSON-type object of results output.
        """
        return FormulaOutput(self.output_results).dump_to_json(print_rounding_order)

    def export_to_json(
        self, filename: str = "default", print_rounding_order: int = 4
    ) -> None:
        """Method to print a final result of calculations in a JSON file.

        Arguments:
            filename (str): filename string (should end with .json)
            print_rounding_order (int): print precision (4 digits by default)
        
        Returns:
            None
        """
        printing = FormulaOutput(self.output_results).export_to_json(
            filename, print_rounding_order
        )
        return
