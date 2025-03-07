from functools import lru_cache

from .formula_parser import ChemicalFormulaParser
from .molar_mass import MolarMassCalculation
from .formula_validator import FormulaValidator
from .chem_output import ChemicalOutput
from .utils import round_dict_content


class ChemicalFormula:
    """A class for operations on a single chemical formula.

    It constructs with a formula string and can calculate
    parsed formula, molar mass, mass percent, atomic percent,
    oxide percent from this string using
    [ChemicalFormulaParser][chemsynthcalc.formula_parser.ChemicalFormulaParser] and
    [MolarMassCalculation][chemsynthcalc.molar_mass.MolarMassCalculation].

    Parameters:
        formula (str): String of chemical formula
        precision (int): Value of rounding precision (8 by default)

    Raise:
        ValueError: if precision <= 0

    Examples:
        >>> ChemicalFormula("H2O")
        H2O
        >>> ChemicalFormula("H2O").molar_mass
        18.015
        >>> ChemicalFormula("H2O").mass_percent
        {'H': 11.19067444, 'O': 88.80932556}
    """

    def __init__(self, formula: str = "", precision: int = 8) -> None:
        if FormulaValidator(formula).validate_formula():
            self.initial_formula: str = formula.replace(" ", "")

        if precision > 0:
            self.precision: int = precision
        else:
            raise ValueError("precision <= 0")

    def __repr__(self) -> str:
        return f"ChemicalFormula({self.formula}, {self.precision})"

    def __str__(self) -> str:
        return self.formula

    @property
    @lru_cache(maxsize=1)
    def formula(self) -> str:
        """
        A string of chemical formula.
        It is made a property to be relatively immutable.

        Returns:
            The formula string

        Examples:
            >>> ChemicalFormula("K2SO4").formula
            K2SO4
        """
        return self.initial_formula

    @property
    @lru_cache(maxsize=1)
    def parsed_formula(self) -> dict[str, float]:
        """
        Formula parsed into dictionary keeping the initial atom order.

        Returns:
            Parsed dictionary representation of formula string created \
            by [ChemicalFormulaParser][chemsynthcalc.formula_parser.ChemicalFormulaParser]

        Examples:
            >>> ChemicalFormula("K2SO4").parsed_formula
            {'K': 2.0, 'S': 1.0, 'O': 4.0}
        """
        parsed: dict[str, float] = ChemicalFormulaParser(self.formula).parse_formula()
        return round_dict_content(parsed, self.precision, plus=3)

    @property
    @lru_cache(maxsize=1)
    def molar_mass(self) -> float:
        """
        Molar mass of the compound.

        Returns:
            The [molar mass](https://en.wikipedia.org/wiki/Molar_mass) \
            of the formula (in g/mol), calculated from \
            parsed the formula using [MolarMassCalculation][chemsynthcalc.molar_mass.MolarMassCalculation] \

        Examples:
            >>> ChemicalFormula("K2SO4").molar_mass
            174.252
        """
        return round(
            MolarMassCalculation(self.parsed_formula).calculate_molar_mass(),
            self.precision,
        )

    @property
    @lru_cache(maxsize=1)
    def mass_percent(self) -> dict[str, float]:
        """
        The percentage of mass of atoms in the formula.

        Returns:
            A mass percent or \
            [relative mass fraction](https://en.wikipedia.org/wiki/Mass_fraction_(chemistry)) \
            of atoms in parsed chemical formula. The values of \
            mass content are in % (with 100% sum), not fraction.

        Examples:
            >>> ChemicalFormula("K2SO4").mass_percent
            {'K': 44.87523816, 'S': 18.39864105, 'O': 36.72612079}
        """
        output = MolarMassCalculation(self.parsed_formula).calculate_mass_percent()
        return round_dict_content(output, self.precision)

    @property
    @lru_cache(maxsize=1)
    def atomic_percent(self) -> dict[str, float]:
        """
        Atomic percents of atoms in the formula.

        Returns:
            An atomic percent or \
            [relative mole fraction](https://en.wikipedia.org/wiki/Mole_fraction) \
            dictionary of atoms in a parsed chemical formula. \
            The values of mole content are in % (with 100% sum), \
            not fraction.

        Examples:
            >>> ChemicalFormula("K2SO4").atomic_percent
            {'K': 28.57142857, 'S': 14.28571429, 'O': 57.14285714}
        """
        output = MolarMassCalculation(self.parsed_formula).calculate_atomic_percent()
        return round_dict_content(output, self.precision)

    @property
    @lru_cache(maxsize=1)
    def oxide_percent(self) -> dict[str, float]:
        """
        Oxide percents of metals in formula.

        Returns:
            An oxide percent or \
            [oxide fraction](https://d32ogoqmya1dw8.cloudfront.net/files/introgeo/studio/examples/minex02.pdf) \
            dictionary of atoms in parsed chemical formula. Oxide types are listed \
            in the [chemsynthcalc.periodic_table][] file and can be changed to \
            any oxide formula. The values of oxide content \
            are in % (with 100% sum), not fraction.

        Examples:
            >>> ChemicalFormula("K2SO4").oxide_percent
            {'K2O': 54.05676836, 'SO3': 45.94323164}
        """
        output = MolarMassCalculation(self.parsed_formula).calculate_oxide_percent()
        return round_dict_content(output, self.precision)

    @property
    @lru_cache(maxsize=1)
    def output_results(self) -> dict[str, object]:
        """
        Dictionary of the calculation result output for class.

        Returns:
            Output dictionary for all properties listed above

        Examples:
            >>> ChemicalFormula("K2SO4").output_results
            {'formula': 'K2SO4', 'parsed formula': {'K': 2.0, 'S': 1.0, 'O': 4.0}, \
            'molar mass': 174.252, 'mass percent': {'K': 44.87523816, 'S': 18.39864105, 'O': 36.72612079}, \
            'atomic percent': {'K': 28.57142857, 'S': 14.28571429, 'O': 57.14285714}, \
            'oxide percent': {'K2O': 54.05676836, 'SO3': 45.94323164}}
        """
        return {
            "formula": self.formula,
            "parsed formula": self.parsed_formula,
            "molar mass": self.molar_mass,
            "mass percent": self.mass_percent,
            "atomic percent": self.atomic_percent,
            "oxide percent": self.oxide_percent,
        }

    def print_results(self, print_precision: int = 4) -> None:
        """
        Print a final result of calculations in stdout.

        Arguments:
            print_precision (int): print precision (4 digits by default)
        """
        ChemicalOutput(
            self.output_results, print_precision, obj=self.__class__.__name__
        ).print_results()

    def to_txt(self, filename: str = "default", print_precision: int = 4) -> None:
        """
        Export the final result of the calculations in a txt file.

        Arguments:
            filename (str): filename string (should end with .txt)
            print_precision (int): print precision (4 digits by default)
        """
        ChemicalOutput(
            self.output_results, print_precision, obj=self.__class__.__name__
        ).write_to_txt(filename)

    def to_json(self, print_precision: int = 4) -> str:
        """
        Serialization of output into JSON object.

        Arguments:
            print_precision (int): print precision (4 digits by default)

        Returns:
            A JSON-type object
        """
        return ChemicalOutput(
            self.output_results, print_precision, obj=self.__class__.__name__
        ).dump_to_json()

    def to_json_file(self, filename: str = "default", print_precision: int = 4) -> None:
        """
        Export a final result of calculations in a JSON file.

        Arguments:
            filename (str): filename string (should end with .json)
            print_precision (int): print precision (4 digits by default)
        """
        ChemicalOutput(
            self.output_results, print_precision, obj=self.__class__.__name__
        ).write_to_json_file(filename)
