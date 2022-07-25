from math import gcd
from functools import reduce

def find_lcm(int_list:int) -> int:
    lcm = 1
    for i in int_list:
        lcm = lcm*i//gcd(lcm, i)
    return lcm

def find_gcd(int_list:int) -> int:
    x = reduce(gcd, int_list)
    return x

def stripe_formula_from_coefficients(formula):
    if not formula[0].isdigit():
        return 1.0, formula
    else:
        coef = []
        for i, symbol in enumerate(formula):
            if symbol.isdigit() or symbol == '.':
                coef.append(symbol)
            else:
                return float("".join(coef)), formula[i:]

def merge_list_of_dicts(list_of_dicts):
    def mergeDict(dict1, dict2):
       ''' Merge dictionaries and keep values of common keys as sum'''
       def merge_two_dicts(dict1, dict2):
            dict3 = dict1.copy()
            dict3.update(dict2)
            return dict3

       dict3 = merge_two_dicts(dict1, dict2)
       for key, value in dict3.items():
           if key in dict1 and key in dict2:
                   dict3[key] = sum([value , dict1[key]])
       return dict3

    final_dict = list_of_dicts[0]
    for dict_i in list_of_dicts[1:]:
        final_dict = mergeDict(final_dict, dict_i)
    return final_dict

def intify_coefficients(coefficients:list) -> list:
    floats = [i for i in coefficients if type(i)==float]
    maximum_lenght = max([len(str(i).split('.')[1]) for i in floats])
    coefficients = [int(i*maximum_lenght*10) for i in coefficients]
    divider = reduce(gcd, coefficients)
    coefficients = [i/divider for i in coefficients]
    return coefficients