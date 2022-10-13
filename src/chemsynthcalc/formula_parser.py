from re import findall, match, compile
from collections import Counter
from .periodic_table import periodic_table
from .chem_errors import (
    NoSuchAtom,
    InvalidCharacter,
    MoreThanOneAdduct,
    BracketsNotPaired,
)

class ChemicalFormulaParser:
    """Parser of chemical formulas.
    
    Methods of this class take string of compound chemical formula
    and turn it into a dict of atoms as keys and their coefficients as values.

    Arguments:
        formula (str): formula string

    Attributes:
        atom_regex (str): regular expression for finding atoms in formula
        coefficient_regex (str): regular expression for atoms amounts in formula
        atom_and_coefficient_regex (str): atom_regex+coefficient_regex
        opener_brackets (str): opener brackets variations
        closer_brackets (str): closer brackets variations
        adduct_symbols (str): symbols for adduct notation (water of crystallization most often)
        allowed_symbols (str): regular expression for all symbols allowed in the formula string
        list_of_atoms (list): list of all atoms in the periodic table
    """

    def __init__(self, formula: str) -> None:
        self.atom_regex: str = "([A-Z][a-z]*)"
        self.coefficient_regex: str = "((\d+(\.\d+)?)*)"
        self.atom_and_coefficient_regex: str = self.atom_regex + self.coefficient_regex
        self.opener_brackets: str = "({["
        self.closer_brackets: str = ")}]"
        self.adduct_symbols: str = "*·•"
        self.allowed_symbols = r"[^A-Za-z0-9.({[)}\]*·•]"
        self.formula: str = formula
        self.list_of_atoms: list = [x[0] for x in periodic_table]

    def __transform_adduct(self, formula: str) -> str:
        """Transform adduct notation in formula to general brackets notation.

        Arguments:
            formula (str): initial formula string
        
        Returns:
            str: New formula with adduct in bracket notation
        """
        transformed_formula = formula
        for i, token in enumerate(formula):
            if token in self.adduct_symbols:
                m = match(self.coefficient_regex, formula[i + 1 :]).group(0)
                if m:
                    coef = str(m)
                else:
                    coef = ""
                transformed_formula = (
                    formula[:i] + "(" + formula[i + 1 + len(coef) :] + ")" + str(coef)
                )
        return transformed_formula

    def __dictify(self, tuples: tuple) -> dict:
        """Transform tuples of tuples to a dict of atoms.

        Arguments:
            tuples (tuple): tuple of tuples of atoms

        Returns:
            dict: dictionary of atoms
        """
        res = dict()
        for atom, n, m, k in tuples:
            try:
                res[atom] += float(n or 1)
            except KeyError:
                res[atom] = float(n or 1)
        return res

    def __fuse(self, mol1: dict, mol2: dict, w: int = 1) -> dict:
        """Fuse 2 dicts representing molecules. 
        
        Arguments:
            mol1 (dict): dict of atoms 1
            mol2 (dict): dict of atoms 2
            w (int): weight

        Returns:
            dict: a new fuse dict
        """
        return {
            atom: (mol1.get(atom, 0) + mol2.get(atom, 0)) * w
            for atom in set(mol1) | set(mol2)
        }

    def __parse(self, formula: str) -> dict:
        """Parse the formula string
       
        Recurse on opening brackets to parse the subpart and
        return on closing ones because it is the end of said subpart.
        Formula is the argument of this method due to the complications
        of self. Constructions in recursive functions.

        Arguments:
            formula (str): formula string

        Returns:
            dict: Return the molecule dict and length of parsed part
        """
        q = []
        mol = {}
        i = 0

        while i < len(formula):
            # Using a classic loop allow for manipulating the cursor
            token = formula[i]

            if token in self.closer_brackets:
                # Check for an index for this part
                m = match(self.coefficient_regex, formula[i + 1 :]).group(0)
                if m != "":
                    weight = float(m)
                    i += len(m)
                else:
                    weight = 1
                submol = self.__dictify(
                    findall(self.atom_and_coefficient_regex, "".join(q))
                )
                return self.__fuse(mol, submol, weight), i

            elif token in self.opener_brackets:
                submol, l = self.__parse(formula[i + 1 :])
                mol = self.__fuse(mol, submol)
                # skip the already read submol
                i += l + 1
            else:
                q.append(token)

            i += 1

        # Fuse in all that's left at base level
        return (
            self.__fuse(
                mol,
                self.__dictify(findall(self.atom_and_coefficient_regex, "".join(q))),
            ),
            i,
        )

    def _is_formula_valid(self) -> list:
        """Checks if the formula string is valid for parsing
        
        If it does not contains any characters that are not
        allowed.

        Returns:
            bool: True if all characters are allowed
            list: List of characters that are not allowed
        
        Examples:
            >>> ChemicalFormulaParser("K2SO4*H2O")._is_formula_valid()
            True
            >>> ChemicalFormulaParser("K2SO4ンH2O")._is_formula_valid()
            ['ン']
        """
        search = compile(self.allowed_symbols).findall
        if search(self.formula):
            return search(self.formula)
        else:
            return True

    def are_brackets_balanced(self) -> bool:
        """Check if all kinds of brackets come in pairs.

        Returns:
            bool: True if all brackets are in pairs
        
        Examples:
            >>> ChemicalFormulaParser("(K2SO4)*H2O").are_brackets_balanced()
            True
            >>> ChemicalFormulaParser("(K2SO4)*H2O").are_brackets_balanced()
            False
        """
        c = Counter(self.formula)
        bracket_counter = c["["] == c["]"] and c["{"] == c["}"] and c["("] == c[")"]

        return bracket_counter

    def is_adduct_one(self) -> bool:
        """ Check if there is only one adduct in formula

        Returns:
            bool: True if there is only one adduct symbol
        
        Examples:
            >>> ChemicalFormulaParser("K2SO4*H2O").is_adduct_one()
            True
            >>> ChemicalFormulaParser("K2SO4*H2O*CO2").is_adduct_one()
            False
        """
        c = Counter(self.formula)
        i = 0
        for adduct in self.adduct_symbols:
            if c[adduct] > 0:
                i += c[adduct]
        if i <= 1:
            return True
        else:
            return False

    def are_atoms_legal(self, parsed:dict) -> None:
        """Checks if all parsed atoms belong to the periodic table.
        
        Arguments:
            parsed (dict): dictionary of parsed atoms
        
        Raises:
            NoSuchAtom: if one of the parsed atoms are not in the periodic table

        Returns:
            None
        
        Examples:
            >>> ChemicalFormulaParser("K2SO4*H2O").are_atoms_legal({'K': 2.0, 'S': 1.0, 'O': 5.0, 'H': 2.0})
            None
            >>> ChemicalFormulaParser("KHuSO4*H2O").are_atoms_legal({'K': 1.0, 'Hu':1.0, 'S': 1.0, 'O': 5.0, 'H': 2.0})
            chemsynthcalc.chem_errors.NoSuchAtom: No atom Hu in the periodic table!
        """
        for atom in list(parsed.keys()):
            if atom not in self.list_of_atoms:
                raise NoSuchAtom("No atom %s in the periodic table!" % atom)
        return

    def parse_formula(self) -> dict:
        """Parse the formula and return a dict with occurrences of each atom.

        Raises:
            InvalidCharacter: if some characters in initial string is not in `allowed_characters`
            MoreThanOneAdduct: if there is more than one adduct symbol (i.e. *)
            BracketsNotPaired: if some bracket type is not coming in pairs
        
        Returns:
            dict: Dictionary of parsed atoms

        Examples:
            >>> ChemicalFormulaParser("H2O").parse_formula()
            {'H': 2.0, 'O': 1.0}
            >>> ChemicalFormulaParser("C2H5OH").parse_formula()
            {'C': 2.0, 'H': 6.0, 'O': 1.0}
            >>> ChemicalFormulaParser("(K0.6Na0.4)2SO4*7H2O").parse_formula()
            {'K': 1.2, 'Na': 0.8, 'S': 1.0, 'O': 11.0, 'H': 14.0}
        """
        if self._is_formula_valid() != True:
            raise InvalidCharacter(
                "Invalid character(s) %s in the formula %s"
                % (self._is_formula_valid(), self.formula)
            )

        if not self.is_adduct_one():
            raise MoreThanOneAdduct("More than one adduct in the formula")

        if not self.are_brackets_balanced():
            raise BracketsNotPaired("The brackers are not balanced ![{]$[&?)]}!]")

        transformed = self.__transform_adduct(self.formula)
        parsed = self.__parse(transformed)[0]

        self.are_atoms_legal(parsed)

        # make an ordered atoms dict
        atoms_list = findall(self.atom_regex, self.formula)
        atoms_dict = dict.fromkeys(atoms_list)
        output = {}
        for atom in atoms_dict.keys():
            output[atom] = parsed.get(atom)
        return output