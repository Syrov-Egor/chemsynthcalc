import numpy as np
from warnings import warn
from re import compile
from functools import cached_property, lru_cache
from .chemical_formula import ChemicalFormula
from .reaction_matrix import ChemicalReactionMatrix
from .reaction_balance import Balancer
from .chem_output import ReactionOutput
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
    try_comb:bool = False
    ) -> None:

        self.rounding_order:int = rounding_order
        self.allowed_symbols:str = r'[^a-zA-Z0-9.({[)}\]*·•=<->→⇄+ ]'
        #separators order is important
        self.possible_reaction_separators:list[str] = ['==', '=', '<->', '->', '<>', '>', '→', '⇄']
        self.reactant_separator:str = '+'
        self.types_of_modes:list[str] = ["force", "check", "balance"]
        self.temp_reaction:str = reaction.replace(" ", "")
        self.algorithm:str = ""
        self.try_comb:bool = try_comb
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
    @lru_cache
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

    @property
    @lru_cache
    def separator(self) -> str:
        '''
        Separator between reactants and products of chemical reaction.
        '''
        if self._is_reaction_string_valid():
            return self._is_reaction_string_valid()
        else:
            raise ValueError("Incorrect reaction")

    @property
    @lru_cache
    def reactants(self) -> list:
        '''
        List of initially sptlitted reactants (left side of the reaction string).
        Formulas are splitted by reactant_separator (+) and includes initial
        coefficients (in case of force or check modes).
        '''
        return self.reaction.split(self.separator)[0].split(self.reactant_separator)

    @property
    @lru_cache
    def products(self) -> list:
        '''
        List of initially sptlitted products (right side of the reaction string).
        Formulas are splitted by reactant_separator (+) and includes initial
        coefficients (in case of force or check modes).
        '''
        return self.reaction.split(self.separator)[1].split(self.reactant_separator)

    @property
    @lru_cache
    def compounds(self) -> list:
        '''
        List of all initially sptlitted products (left side and right side).
        '''
        return self.reactants+self.products

    @property
    @lru_cache
    def initial_coefficients(self) -> list:
        '''
        List of initial coefficients striped from the compounds 
        in case they have been entered (generally the case for 
        force and check modes).
        '''
        return [stripe_formula_from_coefficients(compound)[0] for compound in self.compounds]
    
    @property
    @lru_cache
    def formulas(self) -> list:
        '''
        Decomposition of list of formulas from the reaction string:
        every formula is striped from coefficient and become
        the ChemicalFormula class object.
        '''
        striped = [stripe_formula_from_coefficients(compound)[1] for compound in self.compounds]
        return [ChemicalFormula(formula) for formula in striped]

    @property
    @lru_cache
    def parsed_formulas(self) -> list:
        '''
        List of parsed formulas of ChemicalFormula objects list.
        '''
        return [compound.parsed_formula for compound in self.formulas]

    @property
    @lru_cache
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

    @property
    @lru_cache
    def reactant_matrix(self) -> np.array:
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
    @lru_cache
    def product_matrix(self) -> np.array:
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
    @lru_cache
    def molar_masses(self) -> list:
        '''
        List of molar masses (in g/mol) calculated from parsed 
        ChemicalFormula objects list.
        '''
        return [compound.molar_mass for compound in self.formulas]

    def to_integer(self, coefficients:list) -> list:
        '''
        Cast a float to integer in a list if it is integer. 
        '''
        return [int(i) if i.is_integer() else i for i in coefficients]

    def balance_reaction(self, algorithm:str = "inv", intify:bool = True, try_comb:bool = False, max_comb:int = 1e8) -> list:
        '''
        High-level function call for all balancing algorithms.
        '''
        if self.mode != "balance":
            raise ValueError("Reaction balancing is only available in balance mode")

        avaliable_algorithms:list[str] = ["inv", "gpinv", "ppinv", "comb"]
        if algorithm not in avaliable_algorithms:
            raise ValueError("There is no algorithm %s! Please choose between inv, gpinv, ppinv and comb algorithms" % algorithm)
        
        balance = Balancer(self.reactant_matrix, self.product_matrix, self.rounding_order, intify, try_comb, max_comb)

        if algorithm == "inv":
            coefficients = balance.calculate_coefficients_inv()
            if coefficients:
                self.algorithm = "inverse"
                return coefficients
            else:
                print("Can't equalize this reaction by inverse algorithm")
                return None

        elif algorithm == "gpinv":
            coefficients = balance.calculate_coefficients_gpinv()
            print(coefficients)
            if coefficients:
                self.algorithm = "general pseudoinverse"
                return coefficients
            else:
                print("Can't equalize this reaction by general pseudoinverse algorithm")
                return None

        elif algorithm == "ppinv":
            coefficients = balance.calculate_coefficients_ppinv()
            if coefficients:
                self.algorithm = "partial pseudoinverse"
                return coefficients
            else:
                print("Can't equalize this reaction by partial pseudoinverse algorithm")
                return None

        elif algorithm == "comb":
            coefficients = balance.calculate_coefficients_comb()
            if coefficients:
                self.algorithm = "combinatorial"
                return coefficients
            else:
                print("Can't equalize this reaction by combinatorial algorithm")
                return None
        
    @cached_property
    def coefficients(self) -> list:
        '''
        Returns coefficients of the chemical reaction. There are 3 possible
        modes that method can run: 
        1) force mode is when coefficients are entered by user in reaction 
        string and the calculation and the calculation takes place regardless 
        of reaction balance (it gives warning if reaction is not balanced); 
        2) check mode is basically the force mode but it will raise a
        error if reaction is not balanced and will not calculate masses); 
        3) balance mode uses one of three auto-balancing aglorithms described
        in detail in `Balancer` class.

        In the fist two cases, the coefficients are just stripped from original
        formulas entered by user. In case of balance mode, coefficients are
        calculated.
        '''
        if self.mode == "force":
            if not Balancer.is_reaction_balanced(self.reactant_matrix, self.product_matrix, self.initial_coefficients):
                warn("This reaction is not balanced. Use the output at your own risk")
                self.algorithm = "user"
            return self.to_integer(self.initial_coefficients)

        elif self.mode == "check":
            if Balancer.is_reaction_balanced(self.reactant_matrix, self.product_matrix, self.initial_coefficients):
                self.algorithm = "user"
                return self.to_integer(self.initial_coefficients)
            else:
                raise ValueError("This reaction is not balanced!")

        elif self.mode == "balance":
            try:
                balance = Balancer(self.reactant_matrix, self.product_matrix, self.rounding_order, True, self.try_comb).calculate_coefficients_auto()
                self.algorithm = balance[1]
                return balance[0]
            except Exception:
                return None

    def coefficients_check(self, coefficients:list) -> bool:
        '''
        Checking the coefficients.
        '''
        if coefficients == None:
            raise TypeError("Coefficients are None")
        elif any(x <= 0 for x in coefficients):
            raise ValueError("0 or -x in coefficients")
        elif len(coefficients) != self.matrix.shape[1]:
            raise ValueError("number of coefficients should be equal to %s" % self.matrix.shape[1])
        
        return True

    @property
    @lru_cache
    def normalized_coefficients(self) -> list:
        '''
        List of coefficients normalized on target compound
        (target coefficient = 1).
        '''
        if self.coefficients_check(self.coefficients):
            normalized_coefficients = [coef/self.coefficients[self.target] for coef in self.coefficients]
            return [int(i) if i.is_integer() else round(i, self.rounding_order) for i in normalized_coefficients]
    
    @property
    def is_balanced(self) -> bool:
        '''
        Returns if the reaction is balanced with current coefficients.
        '''
        balance = Balancer.is_reaction_balanced(self.reactant_matrix, self.product_matrix, self.coefficients)
        return balance

    def generate_final_reaction(self, coefs:list) -> str:
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

    @property
    @lru_cache
    def final_reaction(self) -> str:
        '''
        Final representasion of reaction with coefficients.
        '''
        return self.generate_final_reaction(self.coefficients)

    @property
    @lru_cache
    def final_reaction_normalized(self) -> str:
        '''
        Final representasion of reaction with normalized coefficients.
        '''
        return self.generate_final_reaction(self.normalized_coefficients)
    
    @property
    @lru_cache
    def masses(self) -> list:
        '''
        List of masses (in grams) of the of formulas in reaction
        calculated with coefficients obtained by any of 3 methods.
        Calculates masses by calculating amount of substance nu (nu=mass/molar mass).
        Coefficients of reaction are normalized to the target. After nu of target compound is
        calculated, it broadcasted to other compound (with respect to its coefficients).
        '''
        if self.coefficients_check(self.coefficients):
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
        search=compile(self.allowed_symbols).search
        if bool(search(self.temp_reaction)):
            return False 
        for separator in self.possible_reaction_separators:
            if self.temp_reaction.find(separator) != -1 and self.temp_reaction.find(self.reactant_separator) != -1:
                if self.temp_reaction.split(separator)[1] != '':
                    return separator
        return False
    
    @property
    @lru_cache 
    def output_results(self) -> dict:
        '''
        Collection of every output of calculated
        chemical reaction properties.
        '''
        output = {
            "initial reaction:" : self.reaction,
            "reaction matrix:" : self.matrix,
            "mode:" : self.mode,
            "formulas:" : self.formulas,
            "coefficients:" : self.coefficients,
            "normalized coefficients:" : self.normalized_coefficients,
            "algorithm:" : self.algorithm,
            "is balanced:" : self.is_balanced,
            "final reaction:" : self.final_reaction,
            "final reaction normalized:" : self.final_reaction_normalized,
            "molar masses:" : self.molar_masses,
            "target:" : self.formulas[self.target],
            "masses:" : self.masses
        } 
        return output
    
    def print_results(self, print_rounding_order:int = 4) -> None:
        '''
        Method to print a final result of calculations
        in terminal.
        '''
        printing = ReactionOutput(self.output_results).print_results(print_rounding_order)
        return
        
    def export_to_txt(self, filename:str='default', print_rounding_order:int = 4)  -> None:
        '''
        Method to print a final result of calculations
        in txt file.
        '''
        printing = ReactionOutput(self.output_results).export_to_txt(filename, print_rounding_order)
        return
        
    def export_to_json(self, filename:str='default', print_rounding_order:int = 4)  -> None:
        '''
        Method to print a final result of calculations
        in JSON file.
        '''
        printing = ReactionOutput(self.output_results).export_to_json(filename, print_rounding_order)
        return