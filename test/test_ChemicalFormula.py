import unittest
import sys
sys.path.append('./src')
from chemsynthcalc import ChemicalFormula

class TestChemicalFormula(unittest.TestCase):

    def setUp(self) -> None:
        self.test_data = [
        ('H2O', {'H': 2, 'O': 1}, 18.015, {'H': 11.19067, 'O': 88.80933}, {'H': 66.66667, 'O': 33.33333}, {'H2O': 100.0}),
        ('K2SO4', {'K': 2.0, 'S': 1.0, 'O': 4.0}, 174.252, {'K': 44.87524, 'S': 18.39864, 'O': 36.72612}, {'K': 28.57143, 'S': 14.28571, 'O': 57.14286}, {'K2O': 54.05677, 'SO3': 45.94323}),
        ('(NH4)2SO4*H2O', {'N': 2.0, 'H': 10.0, 'S': 1.0, 'O': 5.0}, 150.149, {'N': 18.65747, 'H': 6.71333, 'S': 21.35212, 'O': 53.27708}, {'N': 11.11111, 'H': 55.55556, 'S': 5.55556, 'O': 27.77778}, {'NO2': 35.0993, 'H2O': 34.36115, 'SO3': 30.53955}),
        ('(K0.6Na0.4)2[S]O4', {'K': 1.2, 'Na': 0.8, 'S': 1.0, 'O': 4.0}, 161.3656, {'K': 29.07534, 'Na': 11.39772, 'S': 19.86793, 'O': 39.65901}, {'K': 17.14286, 'Na': 11.42857, 'S': 14.28571, 'O': 57.14286}, {'K2O': 35.02419, 'Na2O': 15.36362, 'SO3': 49.61219})
        ]
    
    # set of tests for formula string validity checks
    def test_invalid_character(self) -> None:
        string = "猫H2O"
        self.assertRaises(ValueError, lambda: ChemicalFormula(string).parsed_formula)

    def test_brackets_balanced(self) -> None:
        string = "(NH42SO4"
        self.assertRaises(ValueError, lambda: ChemicalFormula(string).parsed_formula)
    
    def test_is_adduct_one(self) -> None:
        string = "(NH4)2SO4*H2O*K2SO4"
        self.assertRaises(ValueError, lambda: ChemicalFormula(string).parsed_formula)
        string2 = "(NH4)2SO4*H2O·K2SO4"
        self.assertRaises(ValueError, lambda: ChemicalFormula(string2).parsed_formula)

    def test_are_atoms_legal(self) -> None:
        string = "Ca3(JO4)2"
        self.assertRaises(ValueError, lambda: ChemicalFormula(string).parsed_formula)

    # set of tests to check ChemicalFormula properties
    def test_parser(self) -> None:
        for item in self.test_data:
            self.formula = ChemicalFormula(item[0])
            self.assertEqual(self.formula.parsed_formula, item[1])

    def test_molar_mass(self) -> None:
        for item in self.test_data:
            self.formula = ChemicalFormula(item[0])
            self.assertEqual(self.formula.molar_mass, item[2])

    def test_mass_percent(self) -> None:
        for item in self.test_data:
            self.formula = ChemicalFormula(item[0])
            self.assertEqual(self.formula.mass_percent, item[3])

    def test_atomic_percent(self) -> None:
        for item in self.test_data:
            self.formula = ChemicalFormula(item[0])
            self.assertEqual(self.formula.atomic_percent, item[4])
    
    def test_oxide_percent(self) -> None:
        for item in self.test_data:
            self.formula = ChemicalFormula(item[0])
            self.assertEqual(self.formula.oxide_percent, item[5])


if __name__ == '__main__':
    unittest.main()