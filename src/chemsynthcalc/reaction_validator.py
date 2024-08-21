import re

from .chem_errors import EmptyReaction, InvalidCharacter, NoSeparator
from .reaction import Reaction

class ReactionValidator(Reaction):

    def _check_empty_reaction(self) -> bool:
        return self.reaction == ""
    
    def _invalid_charachers(self) -> list[str]:
        return re.compile(self.allowed_symbols).findall(self.reaction)
    
    @staticmethod
    def extract_separator(reaction: str, possible_reaction_separators: list[str]) -> str:
        for separator in possible_reaction_separators:
            if reaction.find(separator) != -1:
                if reaction.split(separator)[1] != "" and reaction.split(separator)[0] != "":
                    return separator
        return ""
    
    def _no_reaction_separator(self) -> bool:
        return self.extract_separator(self.reaction, self.possible_reaction_separators) == ""
    
    def _no_reactant_separator(self) -> bool:
        return self.reaction.find(self.reactant_separator) == -1

    def validate_reaction(self) -> bool:
        if self._check_empty_reaction():
            raise EmptyReaction
        elif self._invalid_charachers():
            raise InvalidCharacter(
                f"Invalid character(s) {self._invalid_charachers()} in reaction"
            )
        elif self._no_reaction_separator():
            raise NoSeparator(f"No separator between reactants and products: {self.possible_reaction_separators}")
        elif self._no_reactant_separator():
            raise NoSeparator(f"No separators between compounds: {self.reactant_separator}")

        return True