import pytest
from chemsynthcalc.formula_validator import FormulaValidator
from chemsynthcalc.chem_errors import (
    EmptyFormula,
    InvalidCharacter,
    NoSuchAtom,
    MoreThanOneAdduct,
    BracketsNotPaired,
)


def test_empty_formula():
    formula = ""
    with pytest.raises(EmptyFormula):
        FormulaValidator(formula).validate_formula()


def test_one_invalid_characher():
    formula = "猫H2O"
    with pytest.raises(InvalidCharacter):
        FormulaValidator(formula).validate_formula()


def test_invalid_atoms():
    formula = "Hu2O"
    with pytest.raises(NoSuchAtom):
        FormulaValidator(formula).validate_formula()


def test_unbalanced_brackets():
    formula = "K2Mg2(SO4)3)"
    with pytest.raises(BracketsNotPaired):
        FormulaValidator(formula).validate_formula()


def test_more_than_one_adduct():
    formula = "K2Mg2(SO4)3*2H2O*HCl"
    with pytest.raises(MoreThanOneAdduct):
        FormulaValidator(formula).validate_formula()


def test_all_valid():
    formula = "K2Mg2(SO4)3*2H2O"
    assert FormulaValidator(formula).validate_formula() == True
