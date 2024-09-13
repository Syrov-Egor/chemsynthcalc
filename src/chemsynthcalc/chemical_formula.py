from functools import lru_cache

from .formula_parser import ChemicalFormulaParser
from .molar_mass import MolarMassCalculation
from .formula_validator import FormulaValidator
from .chem_output import ChemicalOutput
from .utils import round_dict_content


class ChemicalFormula:
    def __init__(self, formula: str = "", precision: int = 8) -> None:
        if FormulaValidator(formula).validate_formula() == True:
            self.formula: str = formula

        if precision > 0:
            self.precision: int = precision
        else:
            raise ValueError("precision <= 0")

    def __repr__(self) -> str:
        return "chemsynthcalc ChemicalFormula object with formula: " + self.formula

    def __str__(self) -> str:
        return self.formula

    @property
    @lru_cache
    def parsed_formula(self) -> dict[str, float]:
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
