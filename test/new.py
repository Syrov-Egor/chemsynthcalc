import sys

sys.path.append('./src')
import numpy as np
from fractions import Fraction
from chemsynthcalc.chemutils import find_gcd, find_lcm
from chemsynthcalc.reaction_balance import Balancer
import chemsynthcalc
from math import gcd
from functools import reduce

from chemsynthcalc import ChemicalFormula, ChemicalReaction



a = "NH4ClO4+HNO3+HCl=HClO4+NOCl+N2O+N2O3+H2O+Cl2"
b = "NH4ClO4+NOCl+Cl2=HClO4+N2O+N2O3+HNO3+HCl+H2O"
c = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
reaction = ChemicalReaction("Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2")
print(reaction.coefficients)
reaction.export_to_json()