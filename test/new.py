from re import A
import sys

sys.path.append('./src')
import numpy as np
from fractions import Fraction
from chemsynthcalc.chemutils import find_gcd, find_lcm
from chemsynthcalc.reaction_balance import Balancer
import chemsynthcalc
import json
from math import gcd
from functools import reduce

from chemsynthcalc import ChemicalFormula, ChemicalReaction

chemsynthcalc.run_test()