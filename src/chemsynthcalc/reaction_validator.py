import re

from .chem_errors import EmptyReaction, InvalidCharacter, NoSeparator
from .reaction import Reaction


class ReactionValidator(Reaction):
    """
    Methods of this class validate the initial input reaction.
    """

    def _check_empty_reaction(self) -> bool:
        """
        Checks if reaction is an empty string.
        """
        return self.reaction == ""

    def _invalid_charachers(self) -> list[str]:
        """
        Checks if reaction string contains invalid characters.

        Returns:
            List of invalid characters
        """
        return re.compile(self.allowed_symbols).findall(self.reaction)

    def _no_reaction_separator(self) -> bool:
        """
        Checks if reaction string contains a separator.
        """
        return self.extract_separator() == ""

    def _no_reactant_separator(self) -> bool:
        """
        Checks if reaction string contains a "+".
        """
        return self.reaction.find(self.reactant_separator) == -1

    def validate_reaction(self) -> bool:
        """
        Validation of the reaction string.
        Calls the private methods of this class in order.

        Raise:
            [EmptyReaction][chemsynthcalc.chem_errors.EmptyReaction] if reaction string is empty. <br />
            [InvalidCharacter][chemsynthcalc.chem_errors.InvalidCharacter] if some characters are invalid. <br />
            [NoSeparator][chemsynthcalc.chem_errors.NoSeparator] if any of separators (reaction or reactant) are missing.

        Returns:
            True if all the checks are OK
        """

        if self._check_empty_reaction():
            raise EmptyReaction
        elif self._invalid_charachers():
            raise InvalidCharacter(
                f"Invalid character(s) {self._invalid_charachers()} in reaction"
            )
        elif self._no_reaction_separator():
            raise NoSeparator(
                f"No separator between reactants and products: {self.possible_reaction_separators}"
            )
        elif self._no_reactant_separator():
            raise NoSeparator(
                f"No separators between compounds: {self.reactant_separator}"
            )

        return True
