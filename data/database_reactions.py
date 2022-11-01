import csv
import sys
import os
import timeit
import time
sys.path.append('./src')
from chemsynthcalc import ChemicalReaction

in_fnam = "data\\text_mined_reactions.csv"
with open(in_fnam, encoding='utf-8') as in_file:
    data = [row[0] for row in csv.reader(in_file)]

def run_mined_io():
    fname = 'data\\text_mined_output.txt'
    with open(fname, 'w') as file:
        sys.stdout = file
        for i, reaction in enumerate(data):
            chemical_reaction = ChemicalReaction(reaction)
            print(i)
            print(chemical_reaction.reaction)
            print(chemical_reaction.coefficients)
            print("Is balanced:", chemical_reaction.is_balanced)
            print("Algorithm:", chemical_reaction.algorithm)
            print("Masses:", chemical_reaction.masses)
    

def run_mined_without_io():
    for i, reaction in enumerate(data):
            chemical_reaction = ChemicalReaction(reaction)
            i
            chemical_reaction.reaction
            chemical_reaction.coefficients
            chemical_reaction.is_balanced
            chemical_reaction.algorithm
            chemical_reaction.masses

CYCLES = 10
time = timeit.timeit(lambda:run_mined_io(), number=CYCLES)
sys.stdout = sys.__stdout__
print(time)
print(time/CYCLES)