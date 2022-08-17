import re
import sys
from wsgiref.handlers import format_date_time
sys.path.append('./src')

from chemsynthcalc import ChemicalFormula, ChemicalReaction
from chemsynthcalc.reaction_balance import Balancer

reaction = ChemicalReaction("RbNO3+La2O3+Nb2O5+NaNO3=RbLaNaNb3O10+NO2+O2")

reaction.print_results(print_rounding_order = 4)

#print(round(1.255, 1))