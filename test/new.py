import sys

sys.path.append('./src')
import numpy as np
import chemsynthcalc
from chemsynthcalc import ChemicalFormula, ChemicalReaction



reaction = ChemicalReaction("H2+O2=H2O", mode = "balance", try_comb=True) 

#reaction.print_results()

formula = ChemicalFormula("(K0.6Na0.4)2[S]O4")
print(reaction.export_to_txt(filename="1.txt"))