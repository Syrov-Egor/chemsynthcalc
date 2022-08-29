import sys
import time
import re
import numpy as np
from warnings import warn
from functools import cached_property
from .chemical_formula import ChemicalFormula
from .reaction_matrix import ChemicalReactionMatrix
from .reaction_balance import Balancer
from .chemutils import stripe_formula_from_coefficients

class ChemicalReaction():
    '''
    A class that represents a chemical reaction and do operations on it.
    Takes a reaction string as an input, an index of target compound (for calculation of masses), 
    mode string (for calculating coefficients), target compound mass and rounding order.
    '''
    def __init__(self, 
    reaction:str = "", 
    target:int = 0, 
    mode:str = "balance",
    target_mass:float = 1.0, 
    rounding_order:int = 8,
    max_comb:int = 1e8) -> None:

        self.rounding_order:int = rounding_order
        #separator order is important
        self.allowed_characters = '[^a-zA-Z0-9.({[)}]*·•=<->→⇄]'
        self.possible_reaction_separators:list[str] = ['==', '=', '<->', '->', '<>', '>', '→', '⇄']
        self.reactant_separator:str = '+'
        self.types_of_modes:list[str] = ["force", "check", "balance", "combinatorial"]
        self.temp_reaction:str = reaction.replace(" ", "")
        self.algorithm:str = ""
        self.max_comb:int = max_comb
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

    @cached_property
    def reaction(self) -> str:
        '''
        Initial reaction string with validity check.
        '''
        if self.temp_reaction == "":
            raise ValueError("No reaction!")

        if self._is_reaction_string_valid():
            return self.temp_reaction
        else:
            raise ValueError("Incorrect reaction")

    @cached_property
    def separator(self) -> str:
        '''
        Separator between reactants and products of chemical reaction.
        '''
        if self._is_reaction_string_valid():
            return self._is_reaction_string_valid()
        else:
            raise ValueError("Incorrect reaction")

    @cached_property
    def reactants(self) -> list:
        '''
        List of initially sptlitted reactants (left side of the reaction string).
        Formulas are splitted by reactant_separator (+) and includes initial
        coefficients (in case of force or check modes).
        '''
        return self.reaction.split(self.separator)[0].split(self.reactant_separator)

    @cached_property
    def products(self) -> list:
        '''
        List of initially sptlitted products (right side of the reaction string).
        Formulas are splitted by reactant_separator (+) and includes initial
        coefficients (in case of force or check modes).
        '''
        return self.reaction.split(self.separator)[1].split(self.reactant_separator)

    @cached_property
    def compounds(self) -> list:
        '''
        List of all initially sptlitted products (left side and right side).
        '''
        return self.reactants+self.products

    @cached_property
    def initial_coefficients(self) -> list:
        '''
        List of initial coefficients striped from the compounds 
        in case they have been entered (generally the case for 
        force and check modes).
        '''
        return [stripe_formula_from_coefficients(compound)[0] for compound in self.compounds]
    
    @cached_property
    def formulas(self) -> list:
        '''
        Decomposition of list of formulas from the reaction string:
        every formula is striped from coefficient and become
        the ChemicalFormula class object.
        '''
        striped = [stripe_formula_from_coefficients(compound)[1] for compound in self.compounds]
        return [ChemicalFormula(formula) for formula in striped]

    @cached_property
    def parsed_formulas(self) -> list:
        '''
        List of parsed formulas of ChemicalFormula objects list.
        '''
        return [compound.parsed_formula for compound in self.formulas]

    @cached_property
    def matrix(self) -> np.array:
        '''
        The first implementation of reaction matrix method is probably
        belongs to [Blakley](https://doi.org/10.1021/ed059p728). In general,
        a matrix of chemical reaction is composed of coefficients of each
        atom in each compound, giving a 2D array. For example, a matrix of
        reaction KMnO4+HCl=MnCl2+Cl2+H2O+KCl is 
        1) K  [1. 0. 0. 0. 0. 1.]
        2) Mn [1. 0. 1. 0. 0. 0.]
        3) H  [0. 1. 0. 0. 2. 0.]
        4) Cl [0. 1. 2. 2. 0. 1.]
        5) O  [4. 0. 0. 0. 1. 0.]

        The order of rows does not matter for future operations.
        The matrix composes naturally from previously parsed formulas.
        '''
        return ChemicalReactionMatrix(self.parsed_formulas).create_reaction_matrix()

    @cached_property
    def reactant_matrix(self) ->  np.array:
        '''
        Left half of the reaction matrix that consists only of
        reactnats. For example, for reaction 
        KMnO4+HCl=MnCl2+Cl2+H2O+KCl it is:
        1) K  [1. 0.]
        2) Mn [1. 0.]
        3) H  [0. 1.]
        4) Cl [0. 1.]
        5) O  [4. 0.]
        '''
        return self.matrix[:, :len(self.reactants)]

    @cached_property
    def product_matrix(self) ->  np.array:
        '''
        Right half of the reaction matrix that consists only of
        products. For example, for reaction 
        KMnO4+HCl=MnCl2+Cl2+H2O+KCl it is:
        1) K  [0. 0. 0. 1.]
        2) Mn [1. 0. 0. 0.]
        3) H  [0. 0. 2. 0.]
        4) Cl [2. 2. 0. 1.]
        5) O  [0. 0. 1. 0.]
        '''
        return self.matrix[:, len(self.reactants):]
    
    @cached_property
    def molar_masses(self) -> list:
        '''
        List of molar masses (in g/mol) calculated from parsed 
        ChemicalFormula objects list.
        '''
        return [compound.molar_mass for compound in self.formulas]

    @cached_property
    def coefficients(self) -> list:
        '''
        Returns coefficients of the chemical reaction. There are 4 possible
        modes that method can run: 
        1) force mode is when coefficients are entered by user in reaction 
        string and the calculation and the calculation takes place regardless 
        of reaction balance (it gives warning if reaction is not balanced); 
        2) check mode is basically the force mode but it will raise a
        error if reaction is not balanced and will not calculate masses); 
        3) balance mode uses one of three auto-balancing aglorithms described
        in detail in `Balancer` class.
        4) combinatorial mode which solves the Diophantine equation
        by enumerating vectors of coefficients.

        In the fist two cases, the coefficients are just stripped from original
        formulas entered by user. In case of balance mode, coefficients are
        calculated.
        '''
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
                balance = Balancer(self.reactant_matrix, self.product_matrix, self.rounding_order).auto_balance_reaction()
                self.algorithm = balance[1]
                return balance[0]
            except Exception:
                return
        
        elif self.mode == "combinatorial":
            try:
                balance = Balancer(self.reactant_matrix, self.product_matrix, self.rounding_order).calculate_coefficients_combinatorial(max_number_of_iterations=self.max_comb)
                self.algorithm = balance[1]
                return balance[0]
            except Exception:
                print("Can't balance this reaction")
                return
    
    @cached_property
    def normalized_coefficients(self) -> list:
        '''
        List of coefficients normalized on target compound
        (target coefficient = 1)
        '''
        if self.coefficients[self.target] != 1.0:
            normalized_coefficients = [coef/self.coefficients[self.target] for coef in self.coefficients]
        return [int(i) if i.is_integer() else round(i, self.rounding_order) for i in normalized_coefficients]

    def generate_final_reaction(self, coefs) -> str:
        '''
        Final reaction string with connotated formulas and 
        calculated coefficients.
        '''
        final_reaction = []
        for i, compound in enumerate(self.formulas):
            if coefs[i] != 1 or coefs[i] != 1.0:
                final_reaction.append(str(coefs[i])+str(compound))
            else:
                final_reaction.append(str(compound))
        final_reaction = (self.reactant_separator).join(final_reaction)
        final_reaction = final_reaction.replace(self.reactants[-1]+self.reactant_separator, self.reactants[-1]+self.separator)
        return final_reaction

    @cached_property
    def final_reaction(self) -> str:
        return self.generate_final_reaction(self.coefficients)

    @cached_property
    def final_reaction_normalized(self) -> str:
        return self.generate_final_reaction(self.normalized_coefficients)
    
    @cached_property
    def masses(self) -> list:
        '''
        List of masses (in grams) of the of formulas in reaction
        calculated with coefficients obtained by any of 4 methods.
        Calculates masses by calculating amount of substance nu (nu=mass/molar mass).
        Coefficients of reaction are normalized to the target. After nu of target compound is
        calculated, it broadcasted to other compound (with respect to its coefficients).
        '''
        nu = self.target_mass/self.molar_masses[self.target]
        masses = [round(molar*nu*self.normalized_coefficients[i], self.rounding_order) 
        for i, molar in enumerate(self.molar_masses)]

        return masses
        
    def _is_reaction_string_valid(self) -> bool:
        '''
        Naively checks if the reaction string is valid for parsing:
        if it contains one of reactants-products separators (listed
        in possible_reaction_separators attribute) AND a 
        reactant_separator (+).
        '''

        for separator in self.possible_reaction_separators:
            if self.temp_reaction.find(separator) != -1 and self.temp_reaction.find(self.reactant_separator) != -1:
                if self.temp_reaction.split(separator)[1] != '':
                    return separator
        return False
    
    def print_results(self, to_file:bool=False, print_rounding_order:int = 4) -> None:
        '''
        Method to print a final result of calculations.
        By default prints into terminal; if tofile flag
        is set to True, prints into file with int 
        timestamp in the filename.
        '''
        def print_result():
            print("initial reaction:", self.reaction)
            print("reaction matrix:")
            print(self.matrix)
            print("mode:", self.mode)
            print("coefficients:", self.coefficients)
            print("normalized coefficients:", self.normalized_coefficients)
            if self.mode == "balance" or self.mode == "combinatorial":
                print("balanced by", self.algorithm)
            if self.mode == "check":
                print("reaction is balanced")
            print("target:", self.formulas[self.target])
            print("final reaction:", self.final_reaction)
            print("final reaction normalized:", self.final_reaction_normalized)
            for i, formula in enumerate(self.formulas):
                print("%s: M = %s g/mol, m = %s g" % (
                    formula,
                    '%.{0}f'.format(print_rounding_order) % round(self.molar_masses[i], print_rounding_order),
                    '%.{0}f'.format(print_rounding_order) % round(self.masses[i], print_rounding_order)))
                    
        filename:str = "chemsynthcalc_output_%s.txt" % time.time_ns()
        orig_stdout = sys.stdout
        if to_file:
            file = open(filename, 'w')
            sys.stdout = file
            print_result()
            file.close()
            sys.stdout = orig_stdout
        else:
            print_result()
        return