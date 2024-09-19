from chemsynthcalc.utils import round_dict_content, to_integer, find_gcd, find_lcm


def test_round_dict_content():
    test_dict = {
        "N": 18.657466916196576,
        "H": 6.713331424118708,
        "S": 21.35212355726645,
        "O": 53.277078102418265,
    }
    assert round_dict_content(test_dict, 4) == {
        "N": 18.6575,
        "H": 6.7133,
        "S": 21.3521,
        "O": 53.2771,
    }


def test_to_integer():
    assert to_integer([1.0, 2.1, 2.0, 5]) == [1, 2.1, 2, 5]


def test_find_gcd():
    assert find_gcd([30, 40, 80, 60]) == 10


def test_find_lcm():
    assert find_lcm([30, 40, 80, 60]) == 240
