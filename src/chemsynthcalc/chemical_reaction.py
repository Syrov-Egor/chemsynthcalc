from .reaction_validator import ReactionValidator

class ChemicalReaction:
    def __init__(self, reaction: str) -> None:
        if ReactionValidator(reaction).validate_reaction():
            self.reaction = reaction