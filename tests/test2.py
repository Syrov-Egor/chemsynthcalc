from cgi import test
import csv
import ast
import sys
from wsgiref.handlers import format_date_time
sys.path.append('./src')

from chemsynthcalc import ChemicalFormula, ChemicalReaction
from chemsynthcalc.reaction_balance import Balancer

with open('./tests/testing_reactions.csv') as csvfile:
    test_data = list(csv.reader(csvfile))[1:]

reaction_obj = ChemicalReaction('H2O2+KMnO4+H2SO4=O2+MnSO4+K2SO4+H2O')
balance = Balancer(reaction_obj.reactant_matrix, reaction_obj.product_matrix, reaction_obj.rounding_order)
coefs = balance.calculate_coefficients_Thorne()
print(balance.intify_coefficients(coefs, 10000))

'''
test_data = [(i[0], ast.literal_eval(i[1])) for i in test_data]
for react in test_data:
    print(react[0])
    reaction = ChemicalReaction(react[0])
    balance = Balancer(reaction.reactant_matrix, reaction.product_matrix, rounding_order=reaction.rounding_order)
    print(balance.calculate_coefficients_Thorne())
'''