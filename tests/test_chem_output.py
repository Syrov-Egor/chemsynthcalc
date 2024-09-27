import os
import time
import glob

import pytest

from chemsynthcalc.chem_output import ChemicalOutput
from chemsynthcalc.chemical_formula import ChemicalFormula
from chemsynthcalc.chemical_reaction import ChemicalReaction


def cleanup_files() -> None:
    files_to_delete: list[str] = [
        "CSC_formula_*",
        "CSC_reaction_*",
    ]
    for file_type in files_to_delete:
        paths: list[str] = glob.glob(f"./{file_type}")
        for filename in paths:
            os.remove(filename)


formula: str = "[Ru(C10H8N2)3]Cl2*6H2O"
formula_content: list[str] = [
    "formula: [Ru(C10H8N2)3]Cl2*6H2O\n",
    "parsed formula: {'Ru': 1.0, 'C': 30.0, 'H': 36.0, 'N': 6.0, 'Cl': 2.0, 'O': 6.0}\n",
    "molar mass: 748.624\n",
    "mass percent: {'Ru': 13.5008, 'C': 48.1323, 'H': 4.8473, 'N': 11.2262, 'Cl': 9.4707, 'O': 12.8227}\n",
    "atomic percent: {'Ru': 1.2346, 'C': 37.037, 'H': 44.4444, 'N': 7.4074, 'Cl': 2.4691, 'O': 7.4074}\n",
    "oxide percent: {'RuO2': 6.0802, 'CO2': 60.3267, 'H2O': 14.8168, 'NO2': 12.6126, 'ClO2': 6.1638}\n",
]
formula_json_content: str = (
    '{"formula": "[Ru(C10H8N2)3]Cl2*6H2O", "parsed formula": {"Ru": 1.0, "C": 30.0, "H": 36.0, "N": 6.0, "Cl": 2.0, "O": 6.0}, "molar mass": 748.624, "mass percent": {"Ru": 13.5008, "C": 48.1323, "H": 4.8473, "N": 11.2262, "Cl": 9.4707, "O": 12.8227}, "atomic percent": {"Ru": 1.2346, "C": 37.037, "H": 44.4444, "N": 7.4074, "Cl": 2.4691, "O": 7.4074}, "oxide percent": {"RuO2": 6.0802, "CO2": 60.3267, "H2O": 14.8168, "NO2": 12.6126, "ClO2": 6.1638}}'
)
formula_output = ChemicalFormula(formula).output_results


def test_formula_output_wrong_precision() -> None:
    with pytest.raises(ValueError):
        ChemicalOutput(formula_output, print_precision=-2, obj="formula")


def test_formula_output_wrong_obj() -> None:
    with pytest.raises(ValueError):
        ChemicalOutput(formula_output, print_precision=4, obj="firmula")


def test_formula_name_generation() -> None:
    file_type = "txt"
    assert (
        ChemicalOutput(formula_output, print_precision=4, obj="formula")._generate_filename(  # type: ignore
            file_type
        )[
            :11
        ]
        == f"CSC_formula_{formula}_{time.time_ns()}.{file_type}"[:11]
    )


def test_formula_file_name_txt() -> None:
    ChemicalOutput(formula_output, print_precision=4, obj="formula").write_to_txt(
        filename="default"
    )


def test_formula_file_name_json() -> None:
    ChemicalOutput(formula_output, print_precision=4, obj="formula").write_to_json_file(
        filename="default"
    )


def test_formula_print() -> None:
    ChemicalFormula(formula).print_results()


def test_formula_txt_export() -> None:
    filename = "CSC_formula_test.txt"
    ChemicalFormula(formula).to_txt(filename, print_precision=4)
    with open(filename) as f:
        data = f.readlines()
    assert data == formula_content


def test_formula_to_json() -> None:
    assert ChemicalFormula(formula).to_json() == formula_json_content


def test_formula_json_file_export() -> None:
    filename = "CSC_formula_test.json"
    ChemicalFormula(formula).to_json_file(filename)
    with open(filename) as f:
        data = f.read()
    assert data == formula_json_content


reaction: str = "KI+H2SO4=I2+H2S+K2SO4+H2O"

reaction_output = ChemicalReaction(reaction).output_results

reaction_content: list[str] = [
    "initial reaction: KI+H2SO4=I2+H2S+K2SO4+H2O\n",
    "reaction matrix:\n",
    " [[1. 0. 0. 0. 2. 0.]\n",
    " [1. 0. 2. 0. 0. 0.]\n",
    " [0. 2. 0. 2. 0. 2.]\n",
    " [0. 1. 0. 1. 1. 0.]\n",
    " [0. 4. 0. 0. 4. 1.]]\n",
    "mode: balance\n",
    "formulas: ['KI', 'H2SO4', 'I2', 'H2S', 'K2SO4', 'H2O']\n",
    "coefficients: [8, 5, 4, 1, 4, 4]\n",
    "normalized coefficients: [2, 1.25, 1, 0.25, 1, 1]\n",
    "algorithm: inverse\n",
    "is balanced: True\n",
    "final reaction: 8KI+5H2SO4=4I2+H2S+4K2SO4+4H2O\n",
    "final reaction normalized: 2KI+1.25H2SO4=I2+0.25H2S+K2SO4+H2O\n",
    "molar masses: [165.998, 98.072, 253.8, 34.076, 174.252, 18.015]\n",
    "target: I2\n",
    "masses: [1.3081, 0.483, 1.0, 0.0336, 0.6866, 0.071]\n",
    "KI: M = 165.9980 g/mol, m = 1.3081 g\n",
    "H2SO4: M = 98.0720 g/mol, m = 0.4830 g\n",
    "I2: M = 253.8000 g/mol, m = 1.0000 g\n",
    "H2S: M = 34.0760 g/mol, m = 0.0336 g\n",
    "K2SO4: M = 174.2520 g/mol, m = 0.6866 g\n",
    "H2O: M = 18.0150 g/mol, m = 0.0710 g\n",
]

reaction_json_content: str = (
    '{"initial reaction": "KI+H2SO4=I2+H2S+K2SO4+H2O", "reaction matrix": "[[1. 0. 0. 0. 2. 0.]\\n [1. 0. 2. 0. 0. 0.]\\n [0. 2. 0. 2. 0. 2.]\\n [0. 1. 0. 1. 1. 0.]\\n [0. 4. 0. 0. 4. 1.]]", "mode": "balance", "formulas": ["KI", "H2SO4", "I2", "H2S", "K2SO4", "H2O"], "coefficients": [8, 5, 4, 1, 4, 4], "normalized coefficients": [2, 1.25, 1, 0.25, 1, 1], "algorithm": "inverse", "is balanced": true, "final reaction": "8KI+5H2SO4=4I2+H2S+4K2SO4+4H2O", "final reaction normalized": "2KI+1.25H2SO4=I2+0.25H2S+K2SO4+H2O", "molar masses": [165.998, 98.072, 253.8, 34.076, 174.252, 18.015], "target": "I2", "masses": [1.3081, 0.483, 1.0, 0.0336, 0.6866, 0.071]}'
)


def test_reaction_name_generation() -> None:
    file_type = "txt"
    assert (
        ChemicalOutput(reaction_output, print_precision=4, obj="reaction")._generate_filename(  # type: ignore
            file_type
        )[
            :11
        ]
        == f"CSC_reaction_{reaction}_{time.time_ns()}.{file_type}"[:11]
    )


def test_reaction_print() -> None:
    ChemicalReaction(reaction).print_results()


def test_reaction_txt_export() -> None:
    filename = "CSC_reaction_test.txt"
    ChemicalReaction(reaction).to_txt(filename, print_precision=4)
    with open(filename) as f:
        data = f.readlines()
    assert data == reaction_content


def test_reaction_to_json() -> None:
    assert ChemicalReaction(reaction).to_json() == reaction_json_content


def test_reaction_json_file_export() -> None:
    filename = "CSC_reaction_test.json"
    ChemicalReaction(reaction).to_json_file(filename)
    with open(filename) as f:
        data = f.read()
    assert data == reaction_json_content


def test_cleanup() -> None:
    cleanup_files()
