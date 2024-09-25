import pytest

from chemsynthcalc.chemical_reaction import ChemicalReaction

reaction = "KI+H2SO4=I2+H2S+K2SO4+H2O"


def test_wrong_precision() -> None:
    with pytest.raises(ValueError):
        ChemicalReaction(reaction, precision=-2)


def test_wrong_target_mass() -> None:
    with pytest.raises(ValueError):
        ChemicalReaction(reaction, target_mass=0)


def test_repr() -> None:
    assert (
        ChemicalReaction(reaction).__repr__()
        == "chemsynthcalc ChemicalReaction object: KI+H2SO4=I2+H2S+K2SO4+H2O"
    )


def test_str() -> None:
    assert ChemicalReaction(reaction).__str__() == "KI+H2SO4=I2+H2S+K2SO4+H2O"


def test_calculated_target() -> None:
    assert ChemicalReaction(reaction, target=-1)._calculated_target == 1  # type: ignore


def test_wrong_target_high() -> None:
    with pytest.raises(IndexError):
        ChemicalReaction(reaction, target=4)._calculated_target  # type: ignore


def test_wrong_target_low() -> None:
    with pytest.raises(IndexError):
        ChemicalReaction(reaction, target=-3)._calculated_target  # type: ignore


def test_molar_masses() -> None:
    assert ChemicalReaction(reaction).molar_masses == [
        165.998,
        98.072,
        253.8,
        34.076,
        174.252,
        18.015,
    ]


def test_coefficients() -> None:
    assert ChemicalReaction(reaction).coefficients == [8, 5, 4, 1, 4, 4]


def test_normalized_coefficients() -> None:
    assert ChemicalReaction(reaction).normalized_coefficients == [
        2,
        1.25,
        1,
        0.25,
        1,
        1,
    ]


def test_is_balanced() -> None:
    ChemicalReaction(reaction).normalized_coefficients
    assert ChemicalReaction(reaction).is_balanced == True


def test_final_reaction() -> None:
    assert ChemicalReaction(reaction).final_reaction == "8KI+5H2SO4=4I2+H2S+4K2SO4+4H2O"


def test_final_reaction_normalized() -> None:
    assert (
        ChemicalReaction(reaction).final_reaction_normalized
        == "2KI+1.25H2SO4=I2+0.25H2S+K2SO4+H2O"
    )


def test_masses() -> None:
    assert ChemicalReaction(reaction).masses == [
        1.30810087,
        0.48301812,
        1.0,
        0.0335658,
        0.6865721,
        0.07098109,
    ]
