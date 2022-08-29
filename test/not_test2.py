import csv
import ast
import sys
sys.path.append('./src')

from chemsynthcalc import ChemicalFormula, ChemicalReaction

formula = "Ca3(JO4)2"
print(ChemicalFormula(formula).parsed_formula)
print(ChemicalFormula(formula).molar_mass)
print(ChemicalFormula(formula).mass_percent)
print(ChemicalFormula(formula).atomic_percent)
print(ChemicalFormula(formula).oxide_percent)