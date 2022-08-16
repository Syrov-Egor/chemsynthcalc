import sys
import time
import numpy as np
from warnings import warn
from .chemical_formula import ChemicalFormula
from .calculate_mass import MassCalculation
from .reaction_matrix import ChemicalReactionMatrix
from .reaction_balance import Balancer
from .chemutils import stripe_formula_from_coefficients

class ChemicalReaction():
    '''
    A class that represents a chemical reaction and do operations on it.
    Takes a reaction string as an input, an index of target compound (for calculation of masses), 
    mode string (for calculating coefficients), target compound mass and rounding order.
    '''
    def __init__(self, reaction:str, target:int = 0, mode:str = "balance", target_mass:float = 1.0, rounding_order:int = 5) -> None:
        self.rounding_order:int = rounding_order
        self.possible_reaction_separators:list[str] = ['=', '<->', '->', '<>', '>']
        self.reactant_separator:str = '+'
        self.types_of_modes:list[str] = ["force", "check", "balance"]
        self.temp_reaction:str = reaction.replace(" ", "")
        self.algorithm = None
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
        '''
        Initial reaction string with validity check.
        '''
        if self._is_reaction_string_valid():
            return self.temp_reaction
        else:
            raise ValueError("Incorrect reaction")

    @property
    def separator(self) -> str:
        '''
        Separator between reactants and products of chemical reaction.
        '''
        if self._is_reaction_string_valid():
            return self._is_reaction_string_valid()
        else:
            raise ValueError("Incorrect reaction")

    @property
    def reactants(self) -> list:
        '''
        List of initially sptlitted reactants (left side of the reaction string).
        Formulas are splitted by reactant_separator (+) and includes initial
        coefficients (in case of force or check modes).
        '''
        return self.reaction.split(self.separator)[0].split(self.reactant_separator)

    @property
    def products(self) -> list:
        '''
        List of initially sptlitted products (right side of the reaction string).
        Formulas are splitted by reactant_separator (+) and includes initial
        coefficients (in case of force or check modes).
        '''
        return self.reaction.split(self.separator)[1].split(self.reactant_separator)

    @property
    def compounds(self) -> list:
        '''
        List of all initially sptlitted products (left side and right side).
        '''
        return self.reactants+self.products

    @property
    def initial_coefficients(self) -> list:
        '''
        List of initial coefficients striped from the compounds 
        in case they have been entered (generally the case for 
        force and check modes).
        '''
        return [stripe_formula_from_coefficients(compound)[0] for compound in self.compounds]
    
    @property
    def formulas(self) -> list:
        '''
        Decomposition of list of formulas from the reaction string:
        every formula is striped from coefficient and become
        the ChemicalFormula class object.
        '''
        striped = [stripe_formula_from_coefficients(compound)[1] for compound in self.compounds]
        return [ChemicalFormula(formula) for formula in striped]

    @property
    def parsed_formulas(self) -> list:
        '''
        List of parsed formulas of ChemicalFormula objects list.
        '''
        return [compound.parsed_formula for compound in self.formulas]

    @property
    def matrix(self) -> np.array:
        '''
        The first implementation of reaction matrix method is probably
        belongs to Blakley (https://doi.org/10.1021/ed059p728). In general,
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

    @property
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

    @property
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
    
    @property
    def molar_masses(self) -> list:
        '''
        List of molar masses (in g/mol) calculated from parsed 
        ChemicalFormula objects list.
        '''
        return [compound.molar_mass for compound in self.formulas]

    @property
    def coefficients(self) -> list:
        '''
        Returns coefficients of the chemical reaction. There are 3 possible
        modes that class can run: 
        1) force mode is when coefficients are entered by user in reaction 
        string and the calculation and the calculation takes place regardless 
        of reaction balance (it gives warning if reaction is not balanced); 
        2) check mode is basically the force mode but it will raise a
        error if reaction is not balanced and will not calculate masses); 
        3) balance mode uses one of three auto-balancing aglorithms described
        in detail in Balancer class.

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
                balance = Balancer(self.reactant_matrix, self.product_matrix).auto_balance_reaction()
                self.algorithm = balance[1]
                return balance[0]
            except:
                print("Can't equalize this reaction")
                return []
        
    @property
    def final_reaction(self) -> str:
        '''
        Final reaction string with connotated formulas and 
        calculated coefficients.
        '''
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
        '''
        List of masses (in grams) of the of formulas in reaction
        calculated with coefficients obtained by any of 3 methods.
        '''
        return MassCalculation(self.coefficients, self.molar_masses, self.target, self.target_mass, self.rounding_order).calculate_masses()
        
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
    
    def print_results(self, to_file:bool=False) -> None:
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
            if self.mode == "balance":
                print("balanced by", self.algorithm)
            if self.mode == "check":
                print("reaction is balanced")
            print("final_reaction:", self.final_reaction)
            print("target:", self.formulas[self.target])
            for i, formula in enumerate(self.formulas):
                print("%s: M = %s g/mol, m = %s g" %(formula, self.molar_masses[i], self.masses[i]))
       
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