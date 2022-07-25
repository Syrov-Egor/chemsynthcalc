import sys
import time
from warnings import warn
from .chemical_formula import ChemicalFormula
from .calculate_mass import MassCalculation
from .reaction_matrix import ChemicalReactionMatrix
from .reaction_balance import Balancer
from .chemutils import stripe_formula_from_coefficients

class ChemicalReaction():
    def __init__(self, reaction:str, target:int = 0, mode:str = "balance", target_mass:float = 1.0, rounding_order:int = 5) -> None:
        self.rounding_order:int = rounding_order
        self.possible_reaction_separators:list[str] = ['=', '<->', '->', '<>', '>']
        self.reactant_separator:str = '+'
        self.types_of_modes:list[str] = ["force", "check", "balance"]
        self.temp_reaction:str = reaction.replace(" ", "")
        if mode in self.types_of_modes:
            self.mode:str = mode
        else:
             raise ValueError("There is no mode %s! Please choose between force, check or balance modes" % mode)
        if target<len(self.products):
            self.target:int = target + len(self.reactants)
        else:
            raise ValueError("Target should be in range of products number")
        if target_mass>0:
            self.target_mass:float = target_mass
        else:
            raise ValueError("Target mass cannot be 0 or lower")

    def __repr__(self):
        return str(self.reaction)

    def __str__(self):
        return str(self.reaction)

    @property
    def reaction(self) -> str:
        if self._is_reaction_string_valid():
            return self.temp_reaction
        else:
            raise ValueError("Incorrect reaction")

    @property
    def separator(self) -> str:
        if self._is_reaction_string_valid():
            return self._is_reaction_string_valid()
        else:
            raise ValueError("Incorrect reaction")

    @property
    def reactants(self) -> list:
        return self.reaction.split(self.separator)[0].split(self.reactant_separator)

    @property
    def products(self) -> list:
        return self.reaction.split(self.separator)[1].split(self.reactant_separator)

    @property
    def compounds(self) -> list:
        return self.reactants+self.products

    @property
    def initial_coefficients(self) -> list:
        return [stripe_formula_from_coefficients(compound)[0] for compound in self.compounds]
    
    @property
    def formulas(self) -> list:
        striped = [stripe_formula_from_coefficients(compound)[1] for compound in self.compounds]
        return [ChemicalFormula(formula) for formula in striped]

    @property
    def parsed_formulas(self) -> list:
        return [compound.parsed_formula for compound in self.formulas]

    @property
    def matrix(self) -> list:
        return ChemicalReactionMatrix(self.parsed_formulas).create_reaction_matrix()

    @property
    def reactant_matrix(self) -> list:
        return self.matrix[:, :len(self.reactants)]

    @property
    def product_matrix(self) -> list:
        return self.matrix[:, len(self.reactants):]
    
    @property
    def molar_masses(self) -> list:
        return [compound.molar_mass for compound in self.formulas]

    @property
    def coefficients(self) -> list:

        if self.mode == "force":
            if not Balancer.is_reaction_balanced(self.reactant_matrix, self.product_matrix, self.initial_coefficients):
                warn("This reaction is not balanced. Use the output at your own risk")
            return [int(i) if i.is_integer() else i for i in self.initial_coefficients]

        elif self.mode == "check":
            if Balancer.is_reaction_balanced(self.reactant_matrix, self.product_matrix, self.initial_coefficients):
                return [int(i) if i.is_integer() else i for i in self.initial_coefficients]
            else:
                raise ValueError("This reaction is not balanced!")

        elif self.mode == "balance":
            try:
                balance = Balancer(self.reactant_matrix, self.product_matrix).auto_balance_reaction()
                self.algorithm = balance[1]
                return balance[0]
            except:
                print("Can't equalize this reaction")
                return []
        
    @property
    def final_reaction(self) -> str:
        final_reaction = []
        coefs = self.coefficients
        for i, compound in enumerate(self.formulas):
            if coefs[i] != 1 or coefs[i] != 1.0:
                final_reaction.append(str(coefs[i])+str(compound))
            else:
                final_reaction.append(str(compound))
        final_reaction = (self.reactant_separator).join(final_reaction)
        final_reaction = final_reaction.replace(self.reactants[-1]+self.reactant_separator, self.reactants[-1]+self.separator)
        return final_reaction
    
    @property
    def masses(self) -> list:
        return MassCalculation(self.coefficients, self.molar_masses, self.target, self.target_mass).calculate_masses()
        
    def _is_reaction_string_valid(self) -> bool:
        for separator in self.possible_reaction_separators:
            if self.temp_reaction.find(separator) != -1 and self.temp_reaction.find(self.reactant_separator) != -1:
                if self.temp_reaction.split(separator)[1] != '':
                    return separator
        return False
    
    def print_results(self, tofile:bool=False) -> None:

        def print_result():
            print("initial reaction:", self.reaction)
            print("reaction matrix:")
            print(self.matrix)
            print("mode:", self.mode)
            print("coefficients:", self.coefficients)
            if self.mode == "balance":
                print("balanced by", self.algorithm)
            if self.mode == "check":
                print("reaction is balanced")
            print("final_reaction:", self.final_reaction)
            print("target:", self.formulas[self.target])
            for i, formula in enumerate(self.formulas):
                print("%s: M = %s g/mol, m = %s g" %(formula, self.molar_masses[i], self.masses[i]))
       
        filename = "chemsynthcalc_output_%s.txt" % time.time_ns()
        orig_stdout = sys.stdout
        if tofile:
            file = open(filename, 'w')
            sys.stdout = file
            print_result()
            file.close()
            sys.stdout = orig_stdout
        else:
            print_result()
        return