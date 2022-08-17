import numpy as np
from itertools import permutations

def combinatorial(number_of_compounds:int, max_coefficient:int) -> list:
    initial_array = np.ones(number_of_compounds, dtype="int32")
    initial_array = [1]*number_of_compounds
    number_of_iterations = max_coefficient**number_of_compounds
    for idx in range(len(initial_array)):
        for iter in range(1, max_coefficient):
            initial_array[idx] += 1
            print(initial_array)
    return number_of_iterations

i = 3
number_of_compounds = 3
cart_array = (np.arange(1, i+1),)*number_of_compounds
permuted = np.array(np.meshgrid(*cart_array)).T.reshape(-1,number_of_compounds)

def product(*args, repeat=1):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = [tuple(pool) for pool in args] * repeat
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
        print(result)
    #for prod in result:
        #yield tuple(prod)


list(product(range(1,10), repeat=8))