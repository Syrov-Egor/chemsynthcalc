"""
Module that contains custom errors for use in `ChemSynthCalc`
"""


class NoSuchAtom(Exception):
    """
    Found atom(s) that are not in the periodic table.
    """

    pass


class InvalidCharacter(Exception):
    """
    Found some characters that do not belong in the
    chemical formula or reaction.
    """

    pass


class MoreThanOneAdduct(Exception):
    """
    There is more than one adduct (*).
    """

    pass


class BracketsNotPaired(Exception):
    """
    Some brackets do not come in pairs.
    """

    pass


class NoSuchMode(Exception):
    """
    Invalid calculation mode detected.
    """

    pass


class NoSuchAlgorithm(Exception):
    """
    Invalid calculation algorithm detected.
    """

    pass


class NoSeparator(Exception):
    """
    No separator was found in the reaction string.
    """

    pass


class ReactionNotBalanced(Exception):
    """
    This reaction is not balanced.
    """

    pass


class ReactantProductDifference(Exception):
    """
    The elements in reaction are not evenly 
    distributed in reactants and products: 
    some of atoms are only in one part of 
    reaction.
    """

    pass


class BadCoeffiecients(Exception):
    """
    The coefficients are not not compliant
    (they have no physical meaning).
    """

    pass
