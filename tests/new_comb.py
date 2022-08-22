import sys
sys.path.append('./src')
from chemsynthcalc import ChemicalReaction
import numpy as np
import math

def convert_size(size_bytes):
   if size_bytes == 0:
       return "0B"
   size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
   i = int(math.floor(math.log(size_bytes, 1024)))
   p = math.pow(1024, i)
   s = round(size_bytes / p, 2)
   return "%s %s" % (s, size_name[i])

#reaction_obj = ChemicalReaction('H2O2+KNO3+H2SO4=K2SO4+NO+H2O+O2', mode='combinatorial')
#print(reaction_obj.coefficients)

'''
number_of_compounds = 3
max_coef = 700
i = max_coef + 1
cart_array = (np.arange(1, i, dtype='ubyte'), )*number_of_compounds
permuted = np.array(np.meshgrid(*cart_array), dtype='ubyte').T.reshape(-1,number_of_compounds)
print(permuted)
print(convert_size(permuted.nbytes))

# 10^8 = 100000000
# x^7 = 100000000


x = 1e8**(1/3)
print(x)
'''

arr = np.array([[1,2,3], [1,2,3]])
print(arr)
broadcasted = np.broadcast_to(arr, [2,3,5])
print(broadcasted)