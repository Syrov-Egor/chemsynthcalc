from math import gcd
from functools import reduce


def round_dict_content(
    input: dict[str, float], precision: int, plus: int = 0
) -> dict[str, float]:
    return {k: round(v, precision + plus) for k, v in input.items()}


def to_integer(coefficients: list[float | int]) -> list[float | int]:
    return [int(i) if i.is_integer() else i for i in coefficients]


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
