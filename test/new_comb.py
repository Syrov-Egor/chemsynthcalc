import sys
sys.path.append('./src')
from chemsynthcalc import ChemicalReaction
import numpy as np
import math
import timeit
import multiprocessing as mp

def convert_size(size_bytes):
   if size_bytes == 0:
       return "0B"
   size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
   i = int(math.floor(math.log(size_bytes, 1024)))
   p = math.pow(1024, i)
   s = round(size_bytes / p, 2)
   return "%s %s" % (s, size_name[i])

reaction_obj = ChemicalReaction('H2O2+KNO3+H2SO4=K2SO4+NO+H2O+O2', mode='combinatorial', max_comb=1e8)
#print(reaction_obj.coefficients)

print(timeit.timeit(lambda: reaction_obj.coefficients))