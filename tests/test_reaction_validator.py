import pytest

from chemsynthcalc.reaction_validator import ReactionValidator
from chemsynthcalc.chem_errors import EmptyReaction, InvalidCharacter, NoSeparator


def test_empty_reaction():
    reaction = ""
    with pytest.raises(EmptyReaction):
        ReactionValidator(reaction).validate_reaction()


def test_invalid_character():
    reaction = "H2+O2=H2カO"
    with pytest.raises(InvalidCharacter):
        ReactionValidator(reaction).validate_reaction()


def test_no_reaction_separator():
    reaction = "H2+O2 H2O"
    with pytest.raises(NoSeparator):
        ReactionValidator(reaction).validate_reaction()


def test_no_reactant_separator():
    reaction = "H2 O2=H2O"
    with pytest.raises(NoSeparator):
        ReactionValidator(reaction).validate_reaction()
