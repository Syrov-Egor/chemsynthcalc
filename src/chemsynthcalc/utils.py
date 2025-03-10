"""
A module with some useful utilities functions.
"""

from math import gcd
from functools import reduce


def round_dict_content(
    input: dict[str, float], precision: int, plus: int = 0
) -> dict[str, float]:
    """
    Round all values of a dictionary to arbitrary precision
    (with optional surplus value).

    Parameters:
        input (dict[str, float]): An input dict
        precision (int): Precision
        plus (int): An optional surplus

    Returns:
        Rounded dictionary
    """
    return {k: round(v, precision + plus) for k, v in input.items()}


def to_integer(coefficients: list[float | int]) -> list[float | int]:
    """
    Cast a float to int if this float is some x.0 (integer), otherwise
    keep a float.

    Parameters:
        coefficients (list[float | int]): Mixed list of floats and ints

    Returns:
        Mixed list of floats and ints
    """
    return [int(i) if i.is_integer() else i for i in coefficients]


def find_lcm(int_list: list[int]) -> int:
    """
    Find Least Common Multiplyer of list of integers

    Parameters:
        int_list (list[int]): A list of integers

    Returns:
        The LCM
    """
    lcm = 1
    for i in int_list:
        lcm = lcm * i // gcd(lcm, i)
    return lcm


def find_gcd(int_list: list[int]) -> int:
    """
    Find Greatest Common Divisor of list of integers

    Parameters:
        int_list (list[int]): A list of integers

    Returns:
        The GCD
    """
    x = reduce(gcd, int_list)
    return x
