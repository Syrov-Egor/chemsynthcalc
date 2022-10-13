from math import gcd
from functools import reduce

"""
A module with different utility functions, 
which are not belong to main classes
"""


def arguments_type_checking(argument, *args) -> None:
    """Checks for type of function argument.
    
    Argments:
        argument (Any): Argument of the function in question
        types (Any): Types that allowed for this argument
    
    Raises:
        TypeError: if agrument type is not one of specified types
    
    Returns:
        None
    """
    if not isinstance(argument, args):
        raise TypeError(
            "The %s parameter is wrong type %s. Should be %s"
            % (argument, type(argument), args)
        )


def find_lcm(int_list: list) -> int:
    """Finds least common multiple of the list of integers.

    Arguments:
        int_list (list): list of integers
    
    Returns:
        int: Least common multiple
    
    Example:
        >>> chemsynthcalc.chemutils.find_lcm([10,20,30]
        60
    """
    lcm = 1
    for i in int_list:
        lcm = lcm * i // gcd(lcm, i)
    return lcm


def find_gcd(int_list: list) -> int:
    """Finds the greatest common divider of the list of integers.

    Arguments:
        int_list (list): list of integers
    
    Returns:
        int: greatest common divider
    
    Example:
        >>> chemsynthcalc.chemutils.find_gcd([10,20,30]
        10
    """
    x = reduce(gcd, int_list)
    return x


def stripe_formula_from_coefficients(formula) -> tuple:
    """A function to stripe the initially parsed formulas from chemical reaction string.
    
    Arguments:
        formula (str): a formula string

    Returns:
        tuple: (float, str) tuple of (coefficient, reaction)
    
    Example:
        >>> chemsynthcalc.chemutils.stripe_formula_from_coefficients("4.3H3PO4")
        (4.3, 'H3PO4')
    """
    if not formula[0].isdigit():
        return 1.0, formula
    else:
        coef = []
        for i, symbol in enumerate(formula):
            if symbol.isdigit() or symbol == ".":
                coef.append(symbol)
            else:
                return float("".join(coef)), formula[i:]

"""
    def intify_coefficients_pinv(self, coefficients:list) -> list:
        '''
        A function to reduce the coefficients to integers by finding 
        greatset common divider.
        '''
        floats = [i for i in coefficients if type(i)==float]
        maximum_lenght = max([len(str(i).split('.')[1]) for i in floats])
        coefficients = [int(i*maximum_lenght*10) for i in coefficients]
        divider = reduce(gcd, coefficients)
        coefficients = [int(i/divider) for i in coefficients]
        return coefficients

"""
