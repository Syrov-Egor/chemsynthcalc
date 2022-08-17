import csv
import ast
import sys
sys.path.append('./src')
from chemsynthcalc import ChemicalFormula
from chemsynthcalc import ChemicalReaction
from chemsynthcalc.reaction_balance import Balancer
import chemsynthcalc
orig_stdout = sys.stdout
#f = open('out.txt', 'w')
#sys.stdout = f

with open('./tests/testing_reactions.csv') as csvfile:
    test_data = list(csv.reader(csvfile))[1:]

test_data = [(i[0], ast.literal_eval(i[1])) for i in test_data]

j=0
for i, reaction in enumerate(test_data):
    print("%s of %s"% (i+1, len(test_data)))
    print(reaction[0])
    print(reaction[1])
    chemical_reaction = ChemicalReaction(reaction[0])
    print(Balancer.is_reaction_balanced(chemical_reaction.reactant_matrix, chemical_reaction.product_matrix, reaction[1]))
    coefficients = chemical_reaction.coefficients
    print(coefficients)
    print(chemical_reaction.algorithm)
    try:
        assert coefficients == reaction[1]
        print('done')
        j+=1
    except:
        print('failed')

print("done %s of %s" % (j, len(test_data)))

sys.stdout = orig_stdout
#f.close()

#reaction = ChemicalReaction("KMnO4+MnSO4+H2O=MnO2+K2SO4+H2SO4")
#balancer = Balancer(reaction.reactant_matrix, reaction.product_matrix, reaction.rounding_order)
#print(balancer.auto_balance_reaction())
