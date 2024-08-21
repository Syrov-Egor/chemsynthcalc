class Reaction:
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