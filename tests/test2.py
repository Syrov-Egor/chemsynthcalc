import re
import sys
from wsgiref.handlers import format_date_time
sys.path.append('./src')

from chemsynthcalc import ChemicalFormula, ChemicalReaction
from chemsynthcalc.reaction_balance import Balancer

reaction_string = "BaCO3 + Y2(CO3)3 + CuCO3 + O2 → YBa2Cu3O7 + CO2"
reaction = ChemicalReaction(
    reaction = reaction_string,
    target = 0,
    target_mass = 3,
    mode = "balance",
)

reaction.print_results(print_rounding_order=4)

#print(round(1.255, 1))