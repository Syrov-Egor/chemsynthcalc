import re
import numpy as np
from re import compile
from functools import cached_property, lru_cache
from .chemical_formula import ChemicalFormula
from .reaction_matrix import ChemicalReactionMatrix
from .reaction_balance import Balancer
from .chem_output import ReactionOutput
from .chemutils import stripe_formula_from_coefficients, arguments_type_checking
from .chem_errors import (
    BadCoeffiecients,
    NoSuchMode,
    NoSeparator,
    InvalidCharacter,
    ReactantProductDifference,
    ReactionNotBalanced,
    NoSuchAlgorithm,
)


class ChemicalReaction:
    """A class that represents a chemical reaction and do operations on it.

    There are three calculation modes:

    1. The "force" mode is used when a user enters coefficients 
    in the reaction string and wants the masses to be calculated 
    whether the reaction is balanced or not.

    2. "check" mode is the same as force, but with reaction
    balance checks.

    3. "balance" mode  tries to automatically calculate
    coefficients from the reaction string.
    
    Important:
        Unlike other properties of this class, the :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.coefficients`
        property can be set directly.

    Arguments:
        reaction (str): a reaction string
        target (int): index of target compound (0 by default, or first compound in the products)
        mode (str): coefficients calculation mode
        target_mass (float): desired mass of target compound (in grams)
        rounding_order (int): value of rounding precision (8 by default)
        try_comb (bool): flag, which determines whether an attempt will be made to equalize the reaction using the combinatorial method in the "balance" mode

    Attributes:
        allowed_symbols (str): characters that allowed in the reaction string
        possible_reaction_separators (list): list of characters that can be used as reactants-products separator
        reactant_separator (str): character that can be used as compounds separator
        types_of_modes (list): allowed types of calculation modes
        algorithm (str): currently used calculation algorithm
    
    Raises:
        ValueError: if rounding order <=0
        NoSuchMode: if `mode` is not in `types_of_modes`
        ValueError: if `target` is out of bounds of product list
        ValueError: if `target` < 0
        ValueError: if `target_mass` <=0

    Examples:
        >>> ChemicalReaction("H2+O2=H2O")
        H2+O2=H2O
        >>> ChemicalReaction("2H2+O2=2H2O", mode="balance").coefficients
        [2, 1, 2]
        >>> ChemicalReaction("2H2+O2=2H2O", mode="check").masses
        [0.11190674, 0.88809326, 1.0]
    """

    def __init__(
        self,
        reaction: str = "",
        target: int = 0,
        mode: str = "balance",
        target_mass: float = 1.0,
        rounding_order: int = 8,
        try_comb: bool = False,
    ) -> None:
        arguments_type_checking(reaction, str)
        arguments_type_checking(target, int)
        arguments_type_checking(mode, str)
        arguments_type_checking(target_mass, int, float)
        arguments_type_checking(rounding_order, int)
        arguments_type_checking(try_comb, bool)

        if rounding_order > 0:
            self.rounding_order: int = rounding_order
        else:
            raise ValueError("rounding order <= 0")

        self.allowed_symbols: str = r"[^a-zA-Z0-9.({[)}\]*·•=<\->→⇄+ ]"
        # separators order is important
        self.possible_reaction_separators: list[str] = [
            "==",
            "=",
            "<->",
            "->",
            "<>",
            ">",
            "→",
            "⇄",
        ]
        self.reactant_separator: str = "+"
        self.types_of_modes: list[str] = ["force", "check", "balance"]
        self.temp_reaction: str = reaction.replace(" ", "")
        self.algorithm: str = "user"
        self.try_comb: bool = try_comb
        if mode in self.types_of_modes:
            self.mode: str = mode
        else:
            raise NoSuchMode(
                "There is no mode %s! Please choose between force, check or balance modes"
                % mode
            )
        if target >= 0:
            if target < len(self.products):
                self.target: int = target + len(self.reactants)
            else:
                raise ValueError("Target should be in range of products number")
        else:
            raise ValueError("Target < 0")

        if target_mass > 0:
            self.target_mass: float = target_mass
        else:
            raise ValueError("Target mass cannot be 0 or lower")

    def __repr__(self):
        return str(self.reaction)

    def __str__(self):
        return str(self.reaction)

    @property
    @lru_cache
    def reaction(self) -> str:
        """Initial reaction

        Returns:
            str: Initial reaction string with validity check.

        Raises:
            ValueError: if the reaction string is empty
            InvalidCharacter: if some of characters in the string are invalid
        
        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").reaction
            KMnO4+HCl=MnCl2+Cl2+H2O+KCl
        """
        if self.temp_reaction == "":
            raise ValueError("No reaction!")

        if self._is_reaction_string_valid() == True:
            return self.temp_reaction
        else:
            raise InvalidCharacter(
                "Invalid character(s) in the reaction string: %s"
                % self._is_reaction_string_valid()
            )

    @property
    @lru_cache
    def separator(self) -> str:
        """Reactants-product separator.

        Returns:
            str: Separator between reactants and products of chemical reaction.
        
        Raises:
            NoSeparator: if no separator was found in the reaction string

        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").separator
            =
        """
        if self._reaction_contains_separator():
            return self._reaction_contains_separator()
        else:
            raise NoSeparator("No separator in reaction")

    @property
    @lru_cache
    def reactants(self) -> list:
        """List of initially split reactants (left side of the reaction string).
        
        Returns:
            list: Formulas on the left side are split by reactant_separator (+) 
            and include initial coefficients (in case of force or check modes).

        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").reactants
            ['KMnO4', 'HCl']
        """
        return self.reaction.split(self.separator)[0].split(self.reactant_separator)

    @property
    @lru_cache
    def products(self) -> list:
        """List of initially split products (right side of the reaction string).
        
        Returns:
            list: Formulas on the right side that have been split by reactant_separator 
            (+) and include initial coefficients (in case of force or check modes).

        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").products
            ['MnCl2', 'Cl2', 'H2O', 'KCl']
        """
        return self.reaction.split(self.separator)[1].split(self.reactant_separator)

    @property
    @lru_cache
    def compounds(self) -> list:
        """List of all initially split products (left side and right side).

        Returns:
            list: Sum of `reactants` and `products`

        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").compounds
            ['KMnO4', 'HCl', 'MnCl2', 'Cl2', 'H2O', 'KCl']
        """
        return self.reactants + self.products

    @property
    @lru_cache
    def initial_coefficients(self) -> list:
        """Initial coefficients

        The coefficients are not equal to 1 in case if 
        they have been entered (generally the case for 
        force and check modes).

        Returns:
            list: List of initial coefficients striped from the compounds.

        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").initial_coefficients
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        """
        return [
            stripe_formula_from_coefficients(compound)[0] for compound in self.compounds
        ]

    @property
    @lru_cache
    def formulas(self) -> list:
        """Decomposition of a list of formulas from the reaction string.

        Returns:
            list: Every formula is striped from coefficient and become
            the :class:`chemsynthcalc.chemical_formula.ChemicalFormula` object.
        
        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").formulas
            [KMnO4, HCl, MnCl2, Cl2, H2O, KCl]
        """
        striped = [
            stripe_formula_from_coefficients(compound)[1] for compound in self.compounds
        ]
        return [ChemicalFormula(formula) for formula in striped]

    @property
    @lru_cache
    def parsed_formulas(self) -> list:
        """List of parsed formulas of :class:`chemsynthcalc.chemical_formula.ChemicalFormula` objects list.
        
        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").parsed_formulas
            [{'K': 1.0, 'Mn': 1.0, 'O': 4.0}, {'H': 1.0, 'Cl': 1.0}, {'Mn': 1.0, 'Cl': 2.0}, {'Cl': 2.0}, {'H': 2.0, 'O': 1.0}, {'K': 1.0, 'Cl': 1.0}]
        """
        return [compound.parsed_formula for compound in self.formulas]

    @property
    @lru_cache
    def matrix(self) -> np.array:
        """Chemical reaction matrix.

        The first implementation of reaction matrix method is probably
        belongs to `Blakley <https://doi.org/10.1021/ed059p728>`_. In general,
        a chemical reaction matrix is composed of the coefficients of each 
        atom in each compound, giving a 2D array. The matrix composes 
        naturally from previously parsed formulas.

        Returns:
            np.array: 2D array of each atom amount in each formula

        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").matrix
            [[1. 0. 0. 0. 0. 1.]  (K)
             [1. 0. 1. 0. 0. 0.]  (Mn)
             [4. 0. 0. 0. 1. 0.]  (O)
             [0. 1. 0. 0. 2. 0.]  (H)
             [0. 1. 2. 2. 0. 1.]] (Cl)
        """
        return ChemicalReactionMatrix(self.parsed_formulas).create_reaction_matrix()

    @property
    @lru_cache
    def reactant_matrix(self) -> np.array:
        """Left half of the reaction matrix
        
        Returns:
            np.array: 2D array of each atom amount in each reactant

        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").reactant_matrix
            [[1. 0.]  (K)
             [1. 0.]  (Mn)
             [4. 0.]  (O)
             [0. 1.]  (H)
             [0. 1.]] (Cl)
        """
        return self.matrix[:, : len(self.reactants)]

    @property
    @lru_cache
    def product_matrix(self) -> np.array:
        """Right half of the reaction matrix
        
        Returns:
            np.array: 2D array of each atom amount in each product

        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").product_matrix
            [[0. 0. 0. 1.]  (K)
             [1. 0. 0. 0.]  (Mn)
             [0. 0. 1. 0.]  (O)
             [0. 0. 2. 0.]  (H)
             [2. 2. 0. 1.]] (Cl)
        """
        return self.matrix[:, len(self.reactants) :]

    @property
    @lru_cache
    def molar_masses(self) -> list:
        """List of molar masses (in g/mol)
        
        Returns:
            list: List of molar masses of each compound in `compounds`
        
        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").molar_masses
            [158.032, 36.458, 125.838, 70.9, 18.015, 74.548]
        """
        return [compound.molar_mass for compound in self.formulas]

    @cached_property
    def coefficients(self) -> list:
        """Coefficients of the chemical reaction.
        
        There are 3 possible modes that the method can run:

        1) force mode is when coefficients are entered by user in reaction
        string and the calculation and the calculation takes place regardless
        of reaction balance (it gives warning if reaction is not balanced);

        2) check mode is basically the force mode but it will raise an
        error if reaction is not balanced and will not calculate masses);

        3) balance mode uses one of three auto-balancing aglorithms described
        in detail in the :class:`chemsynthcalc.reaction_balance.Balancer` class.

        In the fisrt two cases, the coefficients are just stripped from original
        formulas entered by user. In case of balance mode, coefficients are
        calculated.

        Important:
            Unlike other properties of this class, this
            property can be set directly.

        Raises:
            ReactionNotBalanced: if reaction is not balanced in "check" mode.

        Returns:
            list: list of reaction coefficients
            None: if it cannot balance the reaction

        Examples:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl", mode="balance").coefficients
            [2, 16, 2, 5, 8, 2]
            >>> reaction = ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl", mode="check")
            >>> reaction.coefficients = [2, 16, 2, 5, 8, 2]
            >>> reaction.coefficients
            [2, 16, 2, 5, 8, 2]
            >>> ChemicalReaction("2H2+2O2=H2O", mode="force").coefficients
            [2, 2, 1]
        """
        self.check_elements_count()

        if self.mode == "force":
            if not Balancer.is_reaction_balanced(
                self.reactant_matrix, self.product_matrix, self.initial_coefficients
            ):
                self.algorithm = "user"
            return self.to_integer(self.initial_coefficients)

        elif self.mode == "check":
            if Balancer.is_reaction_balanced(
                self.reactant_matrix, self.product_matrix, self.initial_coefficients
            ):
                self.algorithm = "user"
                return self.to_integer(self.initial_coefficients)
            else:
                raise ReactionNotBalanced("This reaction is not balanced!")

        elif self.mode == "balance":
            try:
                balance = Balancer(
                    self.reactant_matrix,
                    self.product_matrix,
                    self.rounding_order,
                    True,
                    self.try_comb,
                ).calculate_coefficients_auto()
                self.algorithm = balance[1]
                return balance[0]
            except Exception:
                return None

    @property
    @lru_cache
    def normalized_coefficients(self) -> list:
        """List of coefficients normalized on target compound.
        
        target coefficient = 1.0

        Returns:
            list: Normalized coefficients.
        
        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").normalized_coefficients
            [1, 8, 1, 2.5, 4, 1]
        """
        if self.coefficients_check(self.coefficients):
            normalized_coefficients = [
                coef / self.coefficients[self.target] for coef in self.coefficients
            ]
            return [
                int(i) if i.is_integer() else round(i, self.rounding_order)
                for i in normalized_coefficients
            ]

    @property
    def is_balanced(self) -> bool:
        """Is the reaction balanced?

        Returns:
            bool: True if the reaction is balanced with current coefficients
        
        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").is_balanced
            True
        """
        balance = Balancer.is_reaction_balanced(
            self.reactant_matrix, self.product_matrix, self.coefficients
        )
        return balance

    @property
    @lru_cache
    def final_reaction(self) -> str:
        """Final representation of the reaction with coefficients.
        
        Returns:
            str: string of the final reaction
        
        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").final_reaction
            "2KMnO4+16HCl=2MnCl2+5Cl2+8H2O+2KCl"
        """
        return self.generate_final_reaction(self.coefficients)

    @property
    @lru_cache
    def final_reaction_normalized(self) -> str:
        """Final representation of the reaction with normalized coefficients.
        
        Returns:
            str: String of the final normalized reaction
        
        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").final_reaction_normalized
            "KMnO4+8HCl=MnCl2+2.5Cl2+4H2O+KCl"
        """
        return self.generate_final_reaction(self.normalized_coefficients)

    @property
    @lru_cache
    def masses(self) -> list:
        """List of masses of compounds (in grams).
        
        List of masses of the of formulas in reaction
        calculated with coefficients obtained by any of the 3 methods.
        Calculates masses by calculating amount of substance nu (nu=mass/molar mass).
        Coefficients of reaction are normalized to the target. After nu of target compound is
        calculated, it broadcasted to other compound (with respect to their coefficients).

        Returns:
            list: List of masses of compound

        Examples:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").masses
            [1.25583687, 2.31777365, 1.0, 1.40855703, 0.57264101, 0.59241247]
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl", target_mass=2.0).masses
            [2.51167374, 4.63554729, 2.0, 2.81711407, 1.14528203, 1.18482493]
        """
        if self.coefficients_check(self.coefficients):
            nu = self.target_mass / self.molar_masses[self.target]
            masses = [
                round(molar * nu * self.normalized_coefficients[i], self.rounding_order)
                for i, molar in enumerate(self.molar_masses)
            ]
            return masses

    @property
    @lru_cache
    def output_results(self) -> dict:
        """Collection of every output of calculated chemical reaction properties.
        
        Returns:
            dict: All outputs collected in one dictionary.
        
        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").output_results
            {'initial reaction': 'KMnO4+HCl=MnCl2+Cl2+H2O+KCl',
            'reaction matrix': array([[1., 0., 0., 0., 0., 1.],
            [1., 0., 1., 0., 0., 0.],
            [4., 0., 0., 0., 1., 0.],
            [0., 1., 0., 0., 2., 0.],
            [0., 1., 2., 2., 0., 1.]]),
            'mode': 'balance',
            'formulas': [KMnO4, HCl, MnCl2, Cl2, H2O, KCl],
            'coefficients': [2, 16, 2, 5, 8, 2],
            'normalized coefficients': [1, 8, 1, 2.5, 4, 1],
            'algorithm': 'inverse', 
            'is balanced': True, 
            'final reaction': '2KMnO4+16HCl=2MnCl2+5Cl2+8H2O+2KCl', 
            'final reaction normalized': 'KMnO4+8HCl=MnCl2+2.5Cl2+4H2O+KCl', 
            'molar masses': [158.032, 36.458, 125.838, 70.9, 18.015, 74.548], 
            'target': MnCl2, 
            'masses': [1.25583687, 2.31777365, 1.0, 1.40855703, 0.57264101, 0.59241247]}
        """
        output = {
            "initial reaction": self.reaction,
            "reaction matrix": self.matrix,
            "mode": self.mode,
            "formulas": self.formulas,
            "coefficients": self.coefficients,
            "normalized coefficients": self.normalized_coefficients,
            "algorithm": self.algorithm,
            "is balanced": self.is_balanced,
            "final reaction": self.final_reaction,
            "final reaction normalized": self.final_reaction_normalized,
            "molar masses": self.molar_masses,
            "target": self.formulas[self.target],
            "masses": self.masses,
        }
        return output

    def _is_reaction_string_valid(self) -> list:
        """Checks if the reaction string is valid for parsing.
        
        The reaction string is valid if it's not contains any 
        characters that are not allowed.

        Returns:
            bool: True if all characters are allowed
            list: List of characters that are not allowed
        
        Examples:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl")._is_reaction_string_valid()
            True
            >>> ChemicalReaction("KMnO4+HCl=фMnCl2+Cl2+H2O+KCl")._is_reaction_string_valid()
            chemsynthcalc.chem_errors.InvalidCharacter: Invalid character(s) in the reaction string: ['ф']
        """
        search = compile(self.allowed_symbols).findall
        if search(self.temp_reaction):
            return search(self.temp_reaction)
        else:
            return True

    def _reaction_contains_separator(self) -> str:
        """Checks if the reaction string contains one of reactants-products separators:

        Returns:
            str: Separator
            None: if no separator was found
        
        Examples:
            >>> ChemicalReaction("KMnO4+HCl⇄MnCl2+Cl2+H2O+KCl")._reaction_contains_separator()
            ⇄
            >>> ChemicalReaction("KMnO4+HClMnCl2+Cl2+H2O+KCl")._reaction_contains_separator()
            chemsynthcalc.chem_errors.NoSeparator: No separator in reaction
        """
        for separator in self.possible_reaction_separators:
            if self.temp_reaction.find(separator) != -1:
                if self.temp_reaction.split(separator)[1] != "":
                    return separator
        return

    def coefficients_check(self, coefficients: list) -> bool:
        """Checking the coefficients.

        Arguments:
            coefficients (list): list of coefficients

        Raises:
            BadCoeffiecients: If Coefficients are None
            BadCoeffiecients: If 0 or -x in coefficients
            BadCoeffiecients: If number of coefficients is not equal to number of columns in reaction matrix
        
        Returns:
            bool: True, if coefficients are good.
        
        Examples:
            >>> reaction = ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl")
            >>> reaction.coefficients = [2, 16, 2, 5, 8, 2]
            >>> reaction.coefficients_check(reaction.coefficients)
            True
            >>> reaction.coefficients = [2, 16, 2, 5, 8]
            >>> reaction.coefficients_check(reaction.coefficients)
            chemsynthcalc.chem_errors.BadCoeffiecients: number of coefficients should be equal to 6
            >>> reaction.coefficients = None
            >>> reaction.coefficients_check(reaction.coefficients)
            chemsynthcalc.chem_errors.BadCoeffiecients: Coefficients are None
            >>> reaction.coefficients = [2, 16, 2, 5, 8, 0]
            >>> reaction.coefficients_check(reaction.coefficients)
            chemsynthcalc.chem_errors.BadCoeffiecients: 0 or -x in coefficients
        """
        if coefficients == None:
            raise BadCoeffiecients("Coefficients are None")
        elif any(x <= 0 for x in coefficients):
            raise BadCoeffiecients("0 or -x in coefficients")
        elif len(coefficients) != self.matrix.shape[1]:
            raise BadCoeffiecients(
                "number of coefficients should be equal to %s" % self.matrix.shape[1]
            )

        return True

    @lru_cache
    def check_elements_count(self) -> set:
        """Checks if all elements are present in both sides of the reaction.

        Does not work in the `force` mode!

        Returns:
            set: Set of atom differences on the left and right sides of the reaction.

        Examples:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").check_elements_count()
            None
            >>> ChemicalReaction("KMnHoO4+HCl=MnCl2+Cl2+H2O+KCl").check_elements_count()
            chemsynthcalc.chem_errors.ReactantProductDifference: Cannot balance this
            reaction, because element(s) {'Ho'} are only in one part of the reaction
        """
        if self.mode != "force":
            reactants = {
                k: v
                for d in self.parsed_formulas[: len(self.reactants)]
                for k, v in d.items()
            }
            products = {
                k: v
                for d in self.parsed_formulas[len(self.reactants) :]
                for k, v in d.items()
            }
            diff = set(reactants.keys()).symmetric_difference(set(products.keys()))
            if diff:
                raise ReactantProductDifference(
                    "Cannot balance this reaction, because element(s) %s are only in one part of the reaction"
                    % diff
                )
            else:
                return
        else:
            return

    def to_integer(self, coefficients: list) -> list:
        """Cast a float to integer in a list if it is integer.

        Arguments:
            coefficients (list): list of coefficients

        Returns:
            list: List of intified coefficients.

        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").to_integer([1.0, 8.0, 1.0, 2.5, 4.0, 1.0])
            [1, 8, 1, 2.5, 4, 1]
        """
        return [int(i) if i.is_integer() else i for i in coefficients]

    def balance_reaction(
        self, algorithm: str = "inv", intify: bool = True, max_comb: int = 1e8
    ) -> list:
        """High-level function call for all balancing algorithms.
        
        Options availiable are:
        1) `"inv"` - matrix inverse Thorne algorithm
        2) `"gpinv"` - general pseudoinverse Risteski algorithm
        3) `"ppinv"` - partial pseudoinverse Risteski algorithm
        4) `"comb"` - combinatorial algorithm

        Important:
            "comb" algorithm works only for integer coefficients
            less than 127 (see `Balancer.comb_algorithm` for details).

        Arguments:
            algorithm (str) - algortihm choice. "inv" by default.
            intify (bool) - determines whether the coefficients should be integers. True by default.
            max_comb (int) - maximum amount of coefficients vectors hat algorithm goes through. 1e8 by default.

        Returns:
            list: List of obtained coefficients.

        Examples:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").balance_reaction(algorithm="inv", intify=True)
            [2, 16, 2, 5, 8, 2]
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").balance_reaction(algorithm="gpinv", intify=False)
            [0.19607843137254966, 1.568627450980392, 0.1960784313725499, 0.49019607843137436, 0.7843137254901958, 0.19607843137254813]
        """
        arguments_type_checking(algorithm, str)
        arguments_type_checking(intify, bool)
        arguments_type_checking(max_comb, int, float)

        self.check_elements_count()

        if self.mode != "balance":
            raise ValueError("Reaction balancing is only available in balance mode")

        avaliable_algorithms: list[str] = ["inv", "gpinv", "ppinv", "comb"]
        if algorithm not in avaliable_algorithms:
            raise NoSuchAlgorithm(
                "There is no algorithm %s! Please choose between inv, gpinv, ppinv and comb algorithms"
                % algorithm
            )

        if max_comb <= 0:
            raise ValueError("Max combinations <= 0")

        balance = Balancer(
            self.reactant_matrix,
            self.product_matrix,
            self.rounding_order,
            intify,
            self.try_comb,
            max_comb,
        )

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

    def generate_final_reaction(self, coefs: list) -> str:
        """Final reaction string with connotated formulas and calculated coefficients.

        Arguments:
            coefs (list): list of coefficients

        Returns:
            str: String of the final reaction
        
        Example:
            >>> ChemicalReaction("KMnO4+HCl=MnCl2+Cl2+H2O+KCl").generate_final_reaction([2, 16, 2, 5, 8, 2])
            "2KMnO4+16HCl=2MnCl2+5Cl2+8H2O+2KCl"
        """
        final_reaction = []
        for i, compound in enumerate(self.formulas):
            if coefs[i] != 1 or coefs[i] != 1.0:
                final_reaction.append(str(coefs[i]) + str(compound))
            else:
                final_reaction.append(str(compound))
        final_reaction = (self.reactant_separator).join(final_reaction)
        final_reaction = final_reaction.replace(
            self.reactants[-1] + self.reactant_separator,
            self.reactants[-1] + self.separator,
        )
        return final_reaction
    
    def print_results(self, print_rounding_order: int = 4) -> None:
        """Method to print a final result of calculations in terminal.

        Arguments:
            print_rounding_order (int): print precision (4 digits by default)
        
        Returns:
            None
        """
        printing = ReactionOutput(self.output_results).print_results(
            print_rounding_order
        )
        return

    def export_to_txt(
        self, filename: str = "default", print_rounding_order: int = 4
    ) -> None:
        """Method to print a final result of calculations in txt file.

        Arguments:
            filename (str): filename string (should end with .txt)
            print_rounding_order (int): print precision (4 digits by default)
        
        Returns:
            None
        """
        printing = ReactionOutput(self.output_results).export_to_txt(
            filename, print_rounding_order
        )
        return

    def as_json(self, print_rounding_order: int = 4) -> str:
        """Serialization of output into JSON object.

        Argruments:
            print_rounding_order (int): print precision (4 digits by default).
        
        Returns:
            str: JSON object of the output results.
        """
        return ReactionOutput(self.output_results).dump_to_json(print_rounding_order)

    def export_to_json(
        self, filename: str = "default", print_rounding_order: int = 4
    ) -> None:
        """Method to print a final result of calculations in JSON file.

        Arguments:
            filename (str): filename string (should end with .json)
            print_rounding_order (int): print precision (4 digits by default)

        Returns:
            None
        """
        printing = ReactionOutput(self.output_results).export_to_json(
            filename, print_rounding_order
        )
        return