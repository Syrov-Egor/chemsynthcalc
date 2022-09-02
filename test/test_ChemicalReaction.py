import unittest
import csv
import ast
import sys
sys.path.append('./src')
from chemsynthcalc import ChemicalReaction

class TestChemicalReaction(unittest.TestCase):
        def setUp(self) -> None:
            with open('./test/testing_reactions.csv') as csvfile:
                reader = list(csv.reader(csvfile))[1:]
            test_data = [(i[0], ast.literal_eval(i[1])) for i in reader]

        # set of tests for reaction string validity checks
        def test_empty_reaction(self) -> None:
            string = ""
            self.assertRaises(ValueError, lambda: ChemicalReaction(string))

        def test_invalid_character(self) -> None:
            string = "H2+O2=фH2O"
            self.assertRaises(ValueError, lambda: ChemicalReaction(string))

        def test_invalid_string(self) -> None:
            string = "H2O2=H2O"
            self.assertRaises(ValueError, lambda: ChemicalReaction(string))

        def test_wrong_reactant_separator(self) -> None:
            string = "H2&O2=H2O"
            self.assertRaises(ValueError, lambda: ChemicalReaction(string))
        
        def test_wrong_separator(self) -> None:
            string = "H2+O2←H2O"
            self.assertRaises(ValueError, lambda: ChemicalReaction(string))

        # set of tests for ChemicalReaction properties
        def test_initial_reaction(self) -> None:
            return
            #self.assertEqual(self.formula.molar_mass, item[2])

if __name__ == '__main__':
    unittest.main()