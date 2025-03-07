class Reaction:
    """ """

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
        for separator in self.possible_reaction_separators:
            if self.reaction.find(separator) != -1:
                if (
                    self.reaction.split(separator)[1] != ""
                    and self.reaction.split(separator)[0] != ""
                ):
                    return separator
        return ""
