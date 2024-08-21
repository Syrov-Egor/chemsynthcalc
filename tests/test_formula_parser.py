import pytest

from chemsynthcalc.formula_parser import ChemicalFormulaParser

parser_test_data: list[tuple[str, dict[str, float]]] = [
    ("(H2O)", {"H": 2.0, "O": 1.0}),
    ("(K0.6Na0.4)2[S]O4", {"K": 1.2, "Na": 0.8, "S": 1.0, "O": 4.0}),
    ("(NH 4)2 SO4*H2O", {"N": 2.0, "H": 10.0, "S": 1.0, "O": 5.0}),
    ("{K2}2Mg2[(SO4)3Ho]2", {"K": 4.0, "Mg": 2.0, "S": 6.0, "O": 24.0, "Ho": 2.0}),
]


@pytest.mark.parametrize("formula,parsed_formula", parser_test_data)
def test_parser(formula: str, parsed_formula: dict[str, float]):
    assert ChemicalFormulaParser(formula).parse_formula() == parsed_formula
