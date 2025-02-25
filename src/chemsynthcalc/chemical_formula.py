from functools import lru_cache

from .formula_parser import ChemicalFormulaParser
from .molar_mass import MolarMassCalculation
from .formula_validator import FormulaValidator
from .chem_output import ChemicalOutput
from .utils import round_dict_content


class ChemicalFormula:
    """A class for operations on a single chemical formula.

    It constructs with a formula string and calculates properties:
    parsed formula, molar mass, mass percent, atomic percent,
    oxide percent from this string using
    :class:`chemsynthcalc.formula_parser.ChemicalFormulaParser` and
    :class:`chemsynthcalc.molar_mass.MolarMassCalculation`.

    Arguments:
        formula (str): string of chemical formula
        precision (int): value of rounding precision (8 by default)

    Raises:
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
            self.formula: str = formula.replace(" ", "")

        if precision > 0:
            self.precision: int = precision
        else:
            raise ValueError("precision <= 0")

    def __repr__(self) -> str:
        return f"ChemicalFormula({self.formula}, {self.precision})"

    def __str__(self) -> str:
        return self.formula

    @property
    @lru_cache
    def parsed_formula(self) -> dict[str, float]:
        """Formula parsed into dict with the initial atom order.

        Returns:
            dict[str, float]: Parsed dictionary representation of formula string created by :class:`chemsynthcalc.formula_parser.ChemicalFormulaParser`.

        Examples:
            >>> ChemicalFormula("K2SO4").parsed_formula
            {'K': 2.0, 'S': 1.0, 'O': 4.0}
        """
        parsed: dict[str, float] = ChemicalFormulaParser(self.formula).parse_formula()
        return round_dict_content(parsed, self.precision, plus=3)

    @property
    @lru_cache
    def molar_mass(self) -> float:
        return round(
            MolarMassCalculation(self.parsed_formula).calculate_molar_mass(),
            self.precision,
        )

    @property
    @lru_cache
    def mass_percent(self) -> dict[str, float]:
        output = MolarMassCalculation(self.parsed_formula).calculate_mass_percent()
        return round_dict_content(output, self.precision)

    @property
    @lru_cache
    def atomic_percent(self) -> dict[str, float]:
        output = MolarMassCalculation(self.parsed_formula).calculate_atomic_percent()
        return round_dict_content(output, self.precision)

    @property
    @lru_cache
    def oxide_percent(self) -> dict[str, float]:
        output = MolarMassCalculation(self.parsed_formula).calculate_oxide_percent()
        return round_dict_content(output, self.precision)

    @property
    @lru_cache
    def output_results(self) -> dict[str, object]:
        return {
            "formula": self.formula,
            "parsed formula": self.parsed_formula,
            "molar mass": self.molar_mass,
            "mass percent": self.mass_percent,
            "atomic percent": self.atomic_percent,
            "oxide percent": self.oxide_percent,
        }

    def print_results(self, print_precision: int = 4) -> None:
        ChemicalOutput(
            self.output_results, print_precision, obj="formula"
        ).print_results()

    def to_txt(self, filename: str = "default", print_precision: int = 4) -> None:
        ChemicalOutput(
            self.output_results, print_precision, obj="formula"
        ).write_to_txt(filename)

    def to_json(self, print_precision: int = 4) -> str:
        return ChemicalOutput(
            self.output_results, print_precision, obj="formula"
        ).dump_to_json()

    def to_json_file(self, filename: str = "default", print_precision: int = 4) -> None:
        ChemicalOutput(
            self.output_results, print_precision, obj="formula"
        ).write_to_json_file(filename)
