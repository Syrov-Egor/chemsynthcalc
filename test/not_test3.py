import sys
sys.path.append('./src')
import numpy as np
from chemsynthcalc import ChemicalReaction
from chemsynthcalc.reaction_balance import Balancer
from chemsynthcalc.chemutils import intify_coefficients

r1  = "PtNH3BrNO3+CuNH3KNO3+BeCO3=C44H50O14.98+Cu(NO3)2+PtO3+Br1.99NO2+K1.97O+BeO+HNO3"
r2  = "PtNH3BrNO3+(Fe(CN)3)4(Fe(CN)2)3+C44H50O15+Cu(NO3)2+K2BeO2=CuNH3KNO3+BeCO3+K3.97Fe(CN)6+PtO2+Br1.96NO2+HNO3"
r3 = "KAu(CN)2+AgRuAuTe8+Fe2(SO4)3+N2Se4+WO3+Na2CO3+H2CO3+HCl=[Ru(C10H8N2)3]Cl2*6H2O+C4H3AuNa1.96OS7+[WCl4(NSeCl)]2+K3.98Fe(CN)6+Au2O3+TeO3+AgO+NO2"
r4 = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
r5 = "H2O2+KMnO4+H2SO4=O2+MnSO4+K2SO4+H2O"
r6 = "H2O2+KNO3+H2SO4=K2SO4+NO+H2O+O2"
r7 = "KMnO4+H2S+H2SO4=S+MnSO4+K2SO4+H2O"

reaction = ChemicalReaction(r6)
print(reaction.matrix)
balance = Balancer(reaction.reactant_matrix, reaction.product_matrix, reaction.rounding_order)
coefs = balance.calculate_coefficients_Risteski()
#coefs = balance.calculate_coefficients_Thorne()
print(coefs)
#print(balance.intify_coefficients(coefs, 100000))
print(intify_coefficients(coefs))