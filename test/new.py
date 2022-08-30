import sys
sys.path.append('./src')
import numpy as np
from chemsynthcalc import ChemicalFormula, ChemicalReaction
from chemsynthcalc.reaction_balance import Balancer


reaction = ChemicalReaction("H2O2+KNO3+H2SO4=K2SO4+NO+H2O+O2", mode = "balance")
#print(reaction.coefficients)
reaction.coefficients = reaction.balance_reaction(algorithm="Combinatorial")
reaction.print_results()