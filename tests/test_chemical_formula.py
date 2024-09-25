import pytest

from chemsynthcalc.chemical_formula import ChemicalFormula

formula: str = "[Ru(C10H8N2)3]Cl2*6H2O"


def test_wrong_precision() -> None:
    with pytest.raises(ValueError):
        ChemicalFormula(formula=formula, precision=-2)


def test_repr() -> None:
    assert (
        ChemicalFormula(formula=formula).__repr__()
        == "chemsynthcalc ChemicalFormula object with formula: [Ru(C10H8N2)3]Cl2*6H2O"
    )


def test_str() -> None:
    assert ChemicalFormula(formula=formula).__str__() == "[Ru(C10H8N2)3]Cl2*6H2O"
