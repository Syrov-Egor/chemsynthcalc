from re import A
import unittest
import csv
import ast
import sys
import json
import os
import numpy as np

sys.path.append("./src")
from chemsynthcalc import ChemicalReaction
from chemsynthcalc.chem_errors import (
    BadCoeffiecients,
    NoSuchMode,
    NoSeparator,
    InvalidCharacter,
    ReactantProductDifference,
    ReactionNotBalanced,
    NoSuchAlgorithm,
)


class TestChemicalReaction(unittest.TestCase):
    def setUp(self) -> None:
        dir_name = os.path.dirname(os.path.realpath(__file__))
        with open(dir_name+"/testing_reactions.csv") as csvfile:
            reader = list(csv.reader(csvfile))[1:]
        self.reactions_set = [(i[0], ast.literal_eval(i[1])) for i in reader]

        self.test_data = [
            {
                "initial reaction": "H2+O2==H2O",
                "reaction matrix": "[[2. 0. 2.]\n [0. 2. 1.]]",
                "mode": "balance",
                "formulas": ["H2", "O2", "H2O"],
                "coefficients": [2, 1, 2],
                "normalized coefficients": [1, 0.5, 1],
                "algorithm": "inverse",
                "is balanced": True,
                "final reaction": "2H2+O2==2H2O",
                "final reaction normalized": "H2+0.5O2==H2O",
                "molar masses": [2.016, 31.998, 18.015],
                "target": "H2O",
                "masses": [0.11190674, 0.88809326, 1.0],
            },
            {
                "initial reaction": "P+HNO3+H2O<->H3PO4+NO",
                "reaction matrix": "[[1. 0. 0. 1. 0.]\n [0. 1. 2. 3. 0.]\n [0. 1. 0. 0. 1.]\n [0. 3. 1. 4. 1.]]",
                "mode": "balance",
                "formulas": ["P", "HNO3", "H2O", "H3PO4", "NO"],
                "coefficients": [3, 5, 2, 3, 5],
                "normalized coefficients": [1, 1.66666667, 0.66666667, 1, 1.66666667],
                "algorithm": "inverse",
                "is balanced": True,
                "final reaction": "3P+5HNO3+2H2O<->3H3PO4+5NO",
                "final reaction normalized": "P+1.66666667HNO3+0.66666667H2O<->H3PO4+1.66666667NO",
                "molar masses": [30.974, 63.012, 18.015, 97.994, 30.006],
                "target": "H3PO4",
                "masses": [0.31608058, 1.07169827, 0.12255852, 1.0, 0.51033737],
            },
            {
                "initial reaction": "Bi2O3+TiO2→Bi4Ti3O12",
                "reaction matrix": "[[ 2.  0.  4.]\n [ 3.  2. 12.]\n [ 0.  1.  3.]]",
                "mode": "balance",
                "formulas": ["Bi2O3", "TiO2", "Bi4Ti3O12"],
                "coefficients": [2, 3, 1],
                "normalized coefficients": [2, 3, 1],
                "algorithm": "inverse",
                "is balanced": True,
                "final reaction": "2Bi2O3+3TiO2→Bi4Ti3O12",
                "final reaction normalized": "2Bi2O3+3TiO2→Bi4Ti3O12",
                "molar masses": [465.957, 79.865, 1171.509],
                "target": "Bi4Ti3O12",
                "masses": [0.79548172, 0.20451828, 1.0],
            },
            {
                "initial reaction": "CaSO4*2H2O⇄CaSO4*0.5H2O+H2O",
                "reaction matrix": "[[1.  1.  0. ]\n [1.  1.  0. ]\n [6.  4.5 1. ]\n [4.  1.  2. ]]",
                "mode": "balance",
                "formulas": ["CaSO4*2H2O", "CaSO4*0.5H2O", "H2O"],
                "coefficients": [2, 2, 3],
                "normalized coefficients": [1, 1, 1.5],
                "algorithm": "inverse",
                "is balanced": True,
                "final reaction": "2CaSO4*2H2O⇄2CaSO4*0.5H2O+3H2O",
                "final reaction normalized": "CaSO4*2H2O⇄CaSO4*0.5H2O+1.5H2O",
                "molar masses": [172.164, 145.1415, 18.015],
                "target": "CaSO4*0.5H2O",
                "masses": [1.18618038, 1.0, 0.18618038],
            },
        ]

    # set of tests for reaction string validity checks
    def test_empty_reaction(self) -> None:
        string = ""
        self.assertRaises(ValueError, lambda: ChemicalReaction(string))

    def test_invalid_character(self) -> None:
        string = "H2+O2=фH2O"
        self.assertRaises(InvalidCharacter, lambda: ChemicalReaction(string))

    def test_invalid_string(self) -> None:
        string = "H2O2H2O"
        self.assertRaises(NoSeparator, lambda: ChemicalReaction(string))

    def test_wrong_reactant_separator(self) -> None:
        string = "H2&O2=H2O"
        self.assertRaises(InvalidCharacter, lambda: ChemicalReaction(string))

    def test_wrong_separator(self) -> None:
        string = "H2+O2←H2O"
        self.assertRaises(InvalidCharacter, lambda: ChemicalReaction(string))

    def test_unbalanceble_reaction(self) -> None:
        string = "Ho+O2=H2O"
        self.assertRaises(
            ReactantProductDifference, lambda: ChemicalReaction(string).coefficients
        )

    # set of tests for ChemicalReaction with wrong parameters
    def test_reaction_type(self) -> None:
        string = 3
        self.assertRaises(TypeError, lambda: ChemicalReaction(string))

    def test_target_type(self) -> None:
        string = "Al2S3+HNO3<>S+NO2+Al(NO3)3+H2O"
        self.assertRaises(TypeError, lambda: ChemicalReaction(string, target="S"))

    def test_mode_type(self) -> None:
        string = "Al2S3+HNO3<>S+NO2+Al(NO3)3+H2O"
        self.assertRaises(TypeError, lambda: ChemicalReaction(string, mode=3))

    def test_target_mass_type(self) -> None:
        string = "Al2S3+HNO3<>S+NO2+Al(NO3)3+H2O"
        self.assertRaises(
            TypeError, lambda: ChemicalReaction(string, target_mass="3.13")
        )

    def test_rounding_order_type(self) -> None:
        string = "Al2S3+HNO3->S+NO2+Al(NO3)3+H2O"
        self.assertRaises(
            TypeError, lambda: ChemicalReaction(string, rounding_order="8")
        )

    def test_try_comb_type(self) -> None:
        string = "Al2S3+HNO3>S+NO2+Al(NO3)3+H2O"
        self.assertRaises(TypeError, lambda: ChemicalReaction(string, try_comb=1))

    def test_balance_algorithm_type(self) -> None:
        string = "Al2S3+HNO3>S+NO2+Al(NO3)3+H2O"
        self.assertRaises(
            TypeError, lambda: ChemicalReaction(string).balance_reaction(algorithm=3)
        )

    def test_balance_intify_type(self) -> None:
        string = "Al2S3+HNO3>S+NO2+Al(NO3)3+H2O"
        self.assertRaises(
            TypeError, lambda: ChemicalReaction(string).balance_reaction(intify=3)
        )

    def test_balance_max_comb_type(self) -> None:
        string = "Al2S3+HNO3>S+NO2+Al(NO3)3+H2O"
        self.assertRaises(
            TypeError, lambda: ChemicalReaction(string).balance_reaction(max_comb="10")
        )

    def test_balance_max_comb_less_than_zero(self) -> None:
        string = "Al2S3+HNO3>S+NO2+Al(NO3)3+H2O"
        self.assertRaises(
            ValueError, lambda: ChemicalReaction(string).balance_reaction(max_comb=-2)
        )

    def test_target_less_than_zero(self) -> None:
        string = "Al2S3+HNO3>S+NO2+Al(NO3)3+H2O"
        self.assertRaises(ValueError, lambda: ChemicalReaction(string, target=-1))

    def test_target_more_than_number_of_compounds(self) -> None:
        string = "Al2S3+HNO3->S+NO2+Al(NO3)3+H2O"
        self.assertRaises(ValueError, lambda: ChemicalReaction(string, target=4))

    def test_wrong_mode(self) -> None:
        string = "Al2S3+HNO3=S+NO2+Al(NO3)3+H2O"
        self.assertRaises(NoSuchMode, lambda: ChemicalReaction(string, mode="calc"))

    def test_target_mass_zero_or_less(self) -> None:
        string = "Al2S3+HNO3=S+NO2+Al(NO3)3+H2O"
        self.assertRaises(ValueError, lambda: ChemicalReaction(string, target_mass=0))

    def test_target_mass_zero_or_less(self) -> None:
        string = "Al2S3+HNO3=S+NO2+Al(NO3)3+H2O"
        self.assertRaises(ValueError, lambda: ChemicalReaction(string, target_mass=0))

    def test_rounding_order_zero_or_less(self) -> None:
        string = "Al2S3+HNO3=S+NO2+Al(NO3)3+H2O"
        self.assertRaises(
            ValueError, lambda: ChemicalReaction(string, rounding_order=0)
        )

    # set of tests for ChemicalReaction properties
    def test_reaction_matrix(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.assertEqual(
                np.array2string(self.reaction.matrix), item.get("reaction matrix")
            )

    def test_mode(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.assertEqual(self.reaction.mode, item.get("mode"))

    def test_formulas(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.assertEqual(
                [str(formula) for formula in self.reaction.formulas],
                item.get("formulas"),
            )

    def test_coefficients(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.assertEqual(self.reaction.coefficients, item.get("coefficients"))

    def test_normalized_coefficients(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.assertEqual(
                self.reaction.normalized_coefficients,
                item.get("normalized coefficients"),
            )

    def test_algorithm(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.reaction.coefficients
            self.assertEqual(self.reaction.algorithm, item.get("algorithm"))

    def test_is_balanced(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.assertEqual(self.reaction.is_balanced, item.get("is balanced"))

    def test_final_reaction(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.assertEqual(self.reaction.final_reaction, item.get("final reaction"))

    def test_final_reaction_normalized(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.assertEqual(
                self.reaction.final_reaction_normalized,
                item.get("final reaction normalized"),
            )

    def test_molar_masses(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.assertEqual(self.reaction.molar_masses, item.get("molar masses"))

    def test_target(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.assertEqual(
                str(self.reaction.formulas[self.reaction.target]), item.get("target")
            )

    def test_masses(self) -> None:
        for item in self.test_data:
            self.reaction = ChemicalReaction(item.get("initial reaction"))
            self.assertEqual(self.reaction.masses, item.get("masses"))

    # test for different calculation modes
    def test_force_mode(self) -> None:
        string = "H2+O2=2H2O"
        self.reaction = ChemicalReaction(string, mode="force")
        self.assertEqual(self.reaction.masses, [0.05595337, 0.88809326, 1.0])

    def test_check_mode_fail(self) -> None:
        string = "H2+O2=2H2O"
        self.reaction = ChemicalReaction(string, mode="check")
        self.assertRaises(ReactionNotBalanced, lambda: self.reaction.coefficients)

    def test_check_mode(self) -> None:
        string = "2H2+O2=2H2O"
        self.reaction = ChemicalReaction(string, mode="check")
        self.assertEqual(self.reaction.masses, [0.11190674, 0.88809326, 1.0])

    def test_balance_method(self) -> None:
        for reaction in self.reactions_set[:107]:
            chemical_reaction = ChemicalReaction(reaction[0], mode="balance")
            self.assertEqual(chemical_reaction.coefficients, reaction[1])

    def test_feeding_coefs_with_zero(self) -> None:
        string = "H2+O2=H2O"
        self.reaction = ChemicalReaction(string, mode="balance")
        self.reaction.coefficients = [0, 1, 2]
        self.assertRaises(BadCoeffiecients, lambda: self.reaction.masses)

    def test_feeding_negative_coefs(self) -> None:
        string = "H2+O2=H2O"
        self.reaction = ChemicalReaction(string, mode="balance")
        self.reaction.coefficients = [-1, 1, 2]
        self.assertRaises(BadCoeffiecients, lambda: self.reaction.masses)

    def test_feeding_wrong_shape_coefs(self) -> None:
        string = "H2+O2=H2O"
        self.reaction = ChemicalReaction(string, mode="balance")
        self.reaction.coefficients = [1, 2]
        self.assertRaises(BadCoeffiecients, lambda: self.reaction.masses)

    # set of tests for reaction balancing methods
    def test_wrong_algorithm(self) -> None:
        chemical_reaction = ChemicalReaction("H2+O2=H2O")
        self.assertRaises(
            NoSuchAlgorithm, lambda: chemical_reaction.balance_reaction(algorithm="ginv")
            )

    def test_inverse_method(self) -> None:
        for reaction in self.reactions_set[:100]:
            chemical_reaction = ChemicalReaction(reaction[0])
            self.assertEqual(
                chemical_reaction.balance_reaction(algorithm="inv"), reaction[1]
            )

    def test_general_pseudoinverse_method(self) -> None:
        for reaction in self.reactions_set[101:107]:
            chemical_reaction = ChemicalReaction(reaction[0])
            self.assertEqual(
                chemical_reaction.balance_reaction(algorithm="gpinv"), reaction[1]
            )

    def test_partial_pseudoinverse_method(self) -> None:
        for reaction in self.reactions_set[108:110]:
            chemical_reaction = ChemicalReaction(reaction[0])
            coefs = chemical_reaction.balance_reaction(algorithm="ppinv", intify=False)
            coefs = [round(x, 10) for x in coefs]
            self.assertEqual(coefs, reaction[1])

    def test_comb_method(self) -> None:
        for reaction in self.reactions_set[111:113]:
            chemical_reaction = ChemicalReaction(reaction[0])
            coefs = chemical_reaction.balance_reaction(algorithm="comb")
            self.assertEqual(coefs, reaction[1])

    # exports tests
    def cleanup_files(self) -> None:
        listdir = os.listdir()
        for item in listdir:
            if item == "CSC_reaction_test.txt" or item == "CSC_reaction_test.json":
                os.remove(item)
        return

    def test_txt_export(self) -> None:
        string = self.test_data[0].get("initial reaction")
        self.formula = ChemicalReaction(string)
        self.formula.export_to_txt(filename="CSC_reaction_test.txt")

    def test_json_serialization(self) -> None:
        string = self.test_data[0].get("initial reaction")
        self.formula = ChemicalReaction(string)
        dictionary = json.loads(self.formula.as_json())
        output = {
            "initial reaction": "H2+O2==H2O",
            "reaction matrix": "[[2. 0. 2.]\n [0. 2. 1.]]",
            "mode": "balance",
            "formulas": ["H2", "O2", "H2O"],
            "coefficients": [2, 1, 2],
            "normalized coefficients": [1, 0.5, 1],
            "algorithm": "inverse",
            "is balanced": True,
            "final reaction": "2H2+O2==2H2O",
            "final reaction normalized": "H2+0.5O2==H2O",
            "molar masses": [2.016, 31.998, 18.015],
            "target": "H2O",
            "masses": [0.1119, 0.8881, 1.0],
        }
        self.assertEqual(dictionary, output)

    def test_json_export(self) -> None:
        string = self.test_data[0].get("initial reaction")
        self.formula = ChemicalReaction(string)
        self.formula.export_to_json(filename="CSC_reaction_test.json")

    def tearDown(self) -> None:
        self.cleanup_files()
        return super().tearDown()


if __name__ == "__main__":
    unittest.main()
