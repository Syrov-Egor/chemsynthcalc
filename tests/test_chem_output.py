import os

import pytest

from chemsynthcalc.chem_output import ChemicalOutput
from chemsynthcalc.chemical_formula import ChemicalFormula


def cleanup_files() -> None:
    listdir = os.listdir()
    files_to_delete = ["CSC_formula_test.txt", "CSC_formula_test.json"]
    for item in listdir:
        if item in files_to_delete:
            os.remove(item)
    return


formula: str = "[Ru(C10H8N2)3]Cl2*6H2O"
content: list[str] = [
    "formula: [Ru(C10H8N2)3]Cl2*6H2O\n",
    "parsed formula: {'Ru': 1.0, 'C': 30.0, 'H': 36.0, 'N': 6.0, 'Cl': 2.0, 'O': 6.0}\n",
    "molar mass: 748.624\n",
    "mass percent: {'Ru': 13.5008, 'C': 48.1323, 'H': 4.8473, 'N': 11.2262, 'Cl': 9.4707, 'O': 12.8227}\n",
    "atomic percent: {'Ru': 1.2346, 'C': 37.037, 'H': 44.4444, 'N': 7.4074, 'Cl': 2.4691, 'O': 7.4074}\n",
    "oxide percent: {'RuO2': 6.0802, 'CO2': 60.3267, 'H2O': 14.8168, 'NO2': 12.6126, 'ClO2': 6.1638}\n",
]
json_content: str = (
    '{"formula": "[Ru(C10H8N2)3]Cl2*6H2O", "parsed formula": {"Ru": 1.0, "C": 30.0, "H": 36.0, "N": 6.0, "Cl": 2.0, "O": 6.0}, "molar mass": 748.624, "mass percent": {"Ru": 13.5008, "C": 48.1323, "H": 4.8473, "N": 11.2262, "Cl": 9.4707, "O": 12.8227}, "atomic percent": {"Ru": 1.2346, "C": 37.037, "H": 44.4444, "N": 7.4074, "Cl": 2.4691, "O": 7.4074}, "oxide percent": {"RuO2": 6.0802, "CO2": 60.3267, "H2O": 14.8168, "NO2": 12.6126, "ClO2": 6.1638}}'
)
output = ChemicalFormula(formula).output_results


def test_formula_wrong_precision() -> None:
    with pytest.raises(ValueError):
        ChemicalOutput(output, print_precision=-2, obj="formula")


def test_formula_wrong_obj() -> None:
    with pytest.raises(ValueError):
        ChemicalOutput(output, print_precision=4, obj="firmula")


def test_formula_txt_export() -> None:
    filename = "CSC_formula_test.txt"
    ChemicalOutput(output, print_precision=4, obj="formula").write_to_txt(filename)
    with open(filename) as f:
        data = f.readlines()
    assert data == content
    cleanup_files()


def test_formula_to_json() -> None:
    assert (
        ChemicalOutput(output, print_precision=4, obj="formula").dump_to_json()
        == json_content
    )


def test_formula_json_file_export() -> None:
    filename = "CSC_formula_test.json"
    ChemicalOutput(output, print_precision=4, obj="formula").write_to_json_file(
        filename
    )
    with open(filename) as f:
        data = f.read()
    assert data == json_content
    cleanup_files()
