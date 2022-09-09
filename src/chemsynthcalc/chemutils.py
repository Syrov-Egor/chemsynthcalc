from math import gcd
from functools import reduce

"""
A module with different utility functions, 
which are not belong to main classes
"""


def arguments_type_checking(argument, *args) -> None:
    """
    Checks for type of function argument.
    """
    if not isinstance(argument, args):
        raise TypeError(
            "The %s parameter is wrong type %s. Should be %s"
            % (argument, type(argument), args)
        )


def find_lcm(int_list: list) -> int:
    """
    Finds least common multiple of the list of integers.
    """
    lcm = 1
    for i in int_list:
        lcm = lcm * i // gcd(lcm, i)
    return lcm


def find_gcd(int_list: list) -> int:
    """
    Finds greatest common divisor of the list of integers.
    """
    x = reduce(gcd, int_list)
    return x


def stripe_formula_from_coefficients(formula) -> tuple:
    """
    A function to stripe the initially parsed formulas from
    chemical reaction string. For example, 4.3H3PO4 -> (4.3, H3PO4).
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


def merge_list_of_dicts(list_of_dicts) -> dict:
    """
    Merge dictionaries and keep values of common keys as sum.
    """

    def mergeDict(dict1, dict2):
        def merge_two_dicts(dict1, dict2):
            dict3 = dict1.copy()
            dict3.update(dict2)
            return dict3

        dict3 = merge_two_dicts(dict1, dict2)
        for key, value in dict3.items():
            if key in dict1 and key in dict2:
                dict3[key] = sum([value, dict1[key]])
        return dict3

    final_dict = list_of_dicts[0]
    for dict_i in list_of_dicts[1:]:
        final_dict = mergeDict(final_dict, dict_i)
    return final_dict


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
