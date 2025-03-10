class Reaction:
    """
    A base class for
    [ReactionDecomposer][chemsynthcalc.reaction_decomposer.ReactionDecomposer] and
    [ReactionValidator][chemsynthcalc.reaction_validator.ReactionValidator]
    containing regexes and symbols.

    Parameters:
        reaction (str): Reaction string

    Attributes:
        allowed_symbols (str): Regex of all symbols allowed in a reaction string
        possible_reaction_separators (list[str]): List of all allowed reaction separators (left and right part separator)
        reactant_separator (str): Only one possible reactant separator ("+")
    """

    def __init__(self, reaction: str) -> None:
        self.allowed_symbols: str = r"[^a-zA-Z0-9.({[)}\]*·•=<\->→⇄+ ]"
        self.possible_reaction_separators: list[str] = [
            "==",
            "=",
            "<->",
            "->",
            "<>",
            ">",
            "→",
            "⇄",
        ]
        self.reactant_separator: str = "+"

        self.reaction = reaction

    def extract_separator(self) -> str:
        """
        Extract one of possible reaction separator from
        the reaction string.

        Returns:
            Separator string if separator is found, empty string if not
        """
        for separator in self.possible_reaction_separators:
            if self.reaction.find(separator) != -1:
                if (
                    self.reaction.split(separator)[1] != ""
                    and self.reaction.split(separator)[0] != ""
                ):
                    return separator
        return ""
