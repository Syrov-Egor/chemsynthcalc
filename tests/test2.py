from cgi import test
import csv
import ast
import sys
sys.path.append('./src')

from chemsynthcalc import ChemicalFormula, ChemicalReaction
from chemsynthcalc.reaction_balance import Balancer

reaction = ChemicalReaction("KAu(CN)2+AgRuAuTe8+Fe2(SO4)3+N2Se4+WO3+Na2CO3+H2CO3+HCl=[Ru(C10H8N2)3]Cl2*6H2O+C4H3AuNa1.96OS7+[WCl4(NSeCl)]2+K3.98Fe(CN)6+Au2O3+TeO3+AgO+NO2")
print(reaction._is_reaction_string_valid())