"""
ChemSynthCalc
=====
This is a Python package for the calculation of the necessary masses
of substances for chemical synthesis directly from the reaction string.
It includes solutions for all intermidiate steps, including chemical
formula parsing, molar mass calculation and reaction balancing with 
different matrix methods.
"""
from .chemical_formula import ChemicalFormula
from .chemical_reaction import ChemicalReaction