from math import gcd
from functools import reduce


def round_dict_content(
    input: dict[str, float], precision: int, plus: int = 0
) -> dict[str, float]:
    return {k: round(v, precision + plus) for k, v in input.items()}


def to_integer(coefficients: list[float | int]) -> list[float | int]:
    return [int(i) if i.is_integer() else i for i in coefficients]


def find_lcm(int_list: list[int]) -> int:
    lcm = 1
    for i in int_list:
        lcm = lcm * i // gcd(lcm, i)
    return lcm


def find_gcd(int_list: list[int]) -> int:
    x = reduce(gcd, int_list)
    return x
