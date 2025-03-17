import pytest

from chemsynthcalc.chemical_reaction import ChemicalReaction
from chemsynthcalc.coefs import Coefficients
from chemsynthcalc.chem_errors import (
    ReactionNotBalanced,
    NoSuchMode,
    BadCoeffiecients,
    ReactantProductDifference,
)


def test_force_mode():
    reaction = "Cr2(SO4)3+Br2+NaOH=NaBr+Na2CrO4+Na2SO4+H2O"
    reaction_obj = ChemicalReaction(reaction, mode="force")
    assert Coefficients(
        reaction_obj.mode,
        reaction_obj.parsed_formulas,
        reaction_obj.matrix,
        reaction_obj.balancer,
        reaction_obj.decomposed_reaction,
    )._calculate_coefficients()[  # type: ignore
        0
    ] == [
        1,
        1,
        1,
        1,
        1,
        1,
        1,
    ]


def test_check_mode_wrong():
    reaction = "Cr2(SO4)3+Br2+NaOH=NaBr+Na2CrO4+Na2SO4+H2O"
    reaction_obj = ChemicalReaction(reaction, mode="check")
    with pytest.raises(ReactionNotBalanced):
        Coefficients(
            reaction_obj.mode,
            reaction_obj.parsed_formulas,
            reaction_obj.matrix,
            reaction_obj.balancer,
            reaction_obj.decomposed_reaction,
        )._calculate_coefficients()  # type: ignore


def test_check_mode_right():
    reaction = "Cr2(SO4)3+3Br2+16NaOH=6NaBr+2Na2CrO4+3Na2SO4+8H2O"
    reaction_obj = ChemicalReaction(reaction, mode="check")
    assert Coefficients(
        reaction_obj.mode,
        reaction_obj.parsed_formulas,
        reaction_obj.matrix,
        reaction_obj.balancer,
        reaction_obj.decomposed_reaction,
    )._calculate_coefficients()[  # type: ignore
        0
    ] == [
        1,
        3,
        16,
        6,
        2,
        3,
        8,
    ]


def test_wrong_mode():
    reaction = "Cr2(SO4)3+Br2+NaOH=NaBr+Na2CrO4+Na2SO4+H2O"
    reaction_obj = ChemicalReaction(reaction, mode="wrongmode")
    with pytest.raises(NoSuchMode):
        Coefficients(
            reaction_obj.mode,
            reaction_obj.parsed_formulas,
            reaction_obj.matrix,
            reaction_obj.balancer,
            reaction_obj.decomposed_reaction,
        )._calculate_coefficients()  # type: ignore


def test_bad_coefs_0():
    reaction = "Cr2(SO4)3+Br2+NaOH=NaBr+Na2CrO4+Na2SO4+H2O"
    reaction_obj = ChemicalReaction(reaction, mode="balance")
    with pytest.raises(BadCoeffiecients):
        Coefficients(
            reaction_obj.mode,
            reaction_obj.parsed_formulas,
            reaction_obj.matrix,
            reaction_obj.balancer,
            reaction_obj.decomposed_reaction,
        ).coefficients_validation(  # type: ignore
            [
                1,
                3,
                16,
                6,
                0,
                3,
                8,
            ]
        )


def test_bad_coefs_len():
    reaction = "Cr2(SO4)3+Br2+NaOH=NaBr+Na2CrO4+Na2SO4+H2O"
    reaction_obj = ChemicalReaction(reaction, mode="balance")
    with pytest.raises(BadCoeffiecients):
        Coefficients(
            reaction_obj.mode,
            reaction_obj.parsed_formulas,
            reaction_obj.matrix,
            reaction_obj.balancer,
            reaction_obj.decomposed_reaction,
        ).coefficients_validation(  # type: ignore
            [
                1,
                3,
                16,
                6,
                2,
                3,
            ]
        )


def test_element_count_validation_left():
    reaction = "Rb2CO3+La2O3+Nb2O5=RbLaNb2O7"
    reaction_obj = ChemicalReaction(reaction, mode="balance")
    with pytest.raises(ReactantProductDifference):
        reaction_obj.coefficients


def test_element_count_validation_right():
    reaction = "Rb2CO3+La2O3+Nb2O5=RbLaNb2O7+CO2+Nd"
    reaction_obj = ChemicalReaction(reaction, mode="balance")
    with pytest.raises(ReactantProductDifference):
        reaction_obj.coefficients


def test_element_count_validation_both():
    reaction = "Rb2CO3+La2O3+Nb2O5=RbLaNb2O7+Nd"
    reaction_obj = ChemicalReaction(reaction, mode="balance")
    with pytest.raises(ReactantProductDifference):
        reaction_obj.coefficients
