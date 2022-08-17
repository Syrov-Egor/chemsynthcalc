import re
import sys
from wsgiref.handlers import format_date_time
sys.path.append('./src')

from chemsynthcalc import ChemicalFormula, ChemicalReaction
from chemsynthcalc.reaction_balance import Balancer

formula = ChemicalFormula("H2O")
print(formula.molar_mass)
reaction = ChemicalReaction("AsH3+KMnO4+H2SO4=H3AsO4+K2SO4+MnSO4+H2O", mode="combinatorial", max_comb_coefficient=12)
reaction.print_results()


#balancer = Balancer(reaction.reactant_matrix, reaction.product_matrix)
#print(Balancer.is_reaction_balanced(reaction.reactant_matrix, reaction.product_matrix, [5,2,3,5,2,1,8]))
#print(balancer.auto_balance_reaction())
