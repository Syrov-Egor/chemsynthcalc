import imp
import unittest
import os
import sys
import json

sys.path.append("./src")
from chemsynthcalc import ChemicalFormula
from chemsynthcalc.chem_errors import (
    NoSuchAtom,
    InvalidCharacter,
    MoreThanOneAdduct,
    BracketsNotPaired,
)


class TestChemicalFormula(unittest.TestCase):
    def setUp(self) -> None:
        self.test_data = [
            {
                "formula": "H2O",
                "parsed formula": {"H": 2.0, "O": 1.0},
                "molar mass": 18.015,
                "mass percent": {"H": 11.19067444, "O": 88.80932556},
                "atomic percent": {"H": 66.66666667, "O": 33.33333333},
                "oxide percent": {"H2O": 100.0},
            },
            {
                "formula": "K2SO4",
                "parsed formula": {"K": 2.0, "S": 1.0, "O": 4.0},
                "molar mass": 174.252,
                "mass percent": {"K": 44.87523816, "S": 18.39864105, "O": 36.72612079},
                "atomic percent": {
                    "K": 28.57142857,
                    "S": 14.28571429,
                    "O": 57.14285714,
                },
                "oxide percent": {"K2O": 54.05676836, "SO3": 45.94323164},
            },
            {
                "formula": "(NH 4)2 SO4*H2O",
                "parsed formula": {"N": 2.0, "H": 10.0, "S": 1.0, "O": 5.0},
                "molar mass": 150.149,
                "mass percent": {
                    "N": 18.65746692,
                    "H": 6.71333142,
                    "S": 21.35212356,
                    "O": 53.2770781,
                },
                "atomic percent": {
                    "N": 11.11111111,
                    "H": 55.55555556,
                    "S": 5.55555556,
                    "O": 27.77777778,
                },
                "oxide percent": {
                    "NO2": 35.09929733,
                    "H2O": 34.36114777,
                    "SO3": 30.5395549,
                },
            },
            {
                "formula": "(K0.6Na0.4)2[S]O4",
                "parsed formula": {"K": 1.2, "Na": 0.8, "S": 1.0, "O": 4.0},
                "molar mass": 161.3656,
                "mass percent": {
                    "K": 29.07534196,
                    "Na": 11.39772046,
                    "S": 19.86792724,
                    "O": 39.65901035,
                },
                "atomic percent": {
                    "K": 17.14285714,
                    "Na": 11.42857143,
                    "S": 14.28571429,
                    "O": 57.14285714,
                },
                "oxide percent": {
                    "K2O": 35.02419351,
                    "Na2O": 15.36362149,
                    "SO3": 49.612185,
                },
            },
        ]

    # set of tests for formula string validity checks
    def test_empty_formula(self) -> None:
        string = ""
        self.assertRaises(ValueError, lambda: ChemicalFormula(string))

    def test_invalid_character(self) -> None:
        string = "猫H2O"
        self.assertRaises(
            InvalidCharacter, lambda: ChemicalFormula(string).parsed_formula
        )

    def test_brackets_balanced(self) -> None:
        string = "(NH42SO4"
        self.assertRaises(
            BracketsNotPaired, lambda: ChemicalFormula(string).parsed_formula
        )

    def test_is_adduct_one(self) -> None:
        string = "(NH4)2SO4*H2O*K2SO4"
        self.assertRaises(
            MoreThanOneAdduct, lambda: ChemicalFormula(string).parsed_formula
        )
        string2 = "(NH4)2SO4*H2O·K2SO4"
        self.assertRaises(
            MoreThanOneAdduct, lambda: ChemicalFormula(string2).parsed_formula
        )

    def test_are_atoms_legal(self) -> None:
        string = "Ca3(JO4)2"
        self.assertRaises(NoSuchAtom, lambda: ChemicalFormula(string).parsed_formula)

    def test_formula_wrong_type(self) -> None:
        string = 3
        self.assertRaises(TypeError, lambda: ChemicalFormula(string))

    def test_rounding_wrong_type(self) -> None:
        string = "H2O"
        self.assertRaises(TypeError, lambda: ChemicalFormula(string, "3"))

    def test_rounding_order_less_than_zero(self) -> None:
        string = "H2O"
        self.assertRaises(
            ValueError, lambda: ChemicalFormula(string, rounding_order=-1)
        )

    # set of tests to check ChemicalFormula properties
    def test_parser(self) -> None:
        for item in self.test_data:
            self.formula = ChemicalFormula(item.get("formula"))
            self.assertEqual(self.formula.parsed_formula, item.get("parsed formula"))

    def test_molar_mass(self) -> None:
        for item in self.test_data:
            self.formula = ChemicalFormula(item.get("formula"))
            self.assertEqual(self.formula.molar_mass, item.get("molar mass"))

    def test_mass_percent(self) -> None:
        for item in self.test_data:
            self.formula = ChemicalFormula(item.get("formula"))
            self.assertEqual(self.formula.mass_percent, item.get("mass percent"))

    def test_atomic_percent(self) -> None:
        for item in self.test_data:
            self.formula = ChemicalFormula(item.get("formula"))
            self.assertEqual(self.formula.atomic_percent, item.get("atomic percent"))

    def test_oxide_percent(self) -> None:
        for item in self.test_data:
            self.formula = ChemicalFormula(item.get("formula"))
            self.assertEqual(self.formula.oxide_percent, item.get("oxide percent"))

    def cleanup_files(self) -> None:
        listdir = os.listdir()
        for item in listdir:
            if item == "CSC_formula_test.txt" or item == "CSC_formula_test.json":
                os.remove(item)
        return

    # exports tests
    def test_txt_export(self) -> None:
        string = self.test_data[0].get("formula")
        self.formula = ChemicalFormula(string)
        self.formula.export_to_txt(filename="CSC_formula_test.txt")

    def test_json_serialization(self) -> None:
        string = self.test_data[0].get("formula")
        self.formula = ChemicalFormula(string)
        dictionary = json.loads(self.formula.as_json())
        output = {
            "formula": "H2O",
            "parsed formula": {"H": 2.0, "O": 1.0},
            "molar mass": 18.015,
            "mass percent": {"H": 11.1907, "O": 88.8093},
            "atomic percent": {"H": 66.6667, "O": 33.3333},
            "oxide percent": {"H2O": 100.0},
        }
        self.assertEqual(dictionary, output)

    def test_json_export(self) -> None:
        string = self.test_data[0].get("formula")
        self.formula = ChemicalFormula(string)
        self.formula.export_to_json(filename="CSC_formula_test.json")

    def tearDown(self) -> None:
        self.cleanup_files()
        return super().tearDown()


if __name__ == "__main__":
    unittest.main()
