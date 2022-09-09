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
    """
    Parser of chemical formulas. Methods of this class take string of compound chemical formula
    and turn it to dict of atoms as keys and their coefficients as values.

    Parameters:
    * `formula:str` - formula string

    For example:
    ```
    >>>ChemicalFormulaParser("H2O").parse_formula()
    {'H': 2.0, 'O': 1.0}
    >>>ChemicalFormulaParser("C2H5OH").parse_formula()
    {'C': 2.0, 'H': 6.0, 'O': 1.0}
    >>>ChemicalFormulaParser("(K0.6Na0.4)2SO4*7H2O").parse_formula()
    {'K': 1.2, 'Na': 0.8, 'S': 1.0, 'O': 11.0, 'H': 14.0}
    ```
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
        """
        Transform adduct notation in formula to general
        brackets notation.
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
        """
        Transform tuples of tuples to a dict of atoms.
        """
        res = dict()
        for atom, n, m, k in tuples:
            try:
                res[atom] += float(n or 1)
            except KeyError:
                res[atom] = float(n or 1)
        return res

    def __fuse(self, mol1: dict, mol2: dict, w: int = 1) -> dict:
        """
        Fuse 2 dicts representing molecules. Return a new dict.
        """
        return {
            atom: (mol1.get(atom, 0) + mol2.get(atom, 0)) * w
            for atom in set(mol1) | set(mol2)
        }

    def __parse(self, formula: str) -> dict:
        """
        Return the molecule dict and length of parsed part.
        Recurse on opening brackets to parse the subpart and
        return on closing ones because it is the end of said subpart.
        Formula is the argument of this method due to the complications
        of self. Constructions in recursive functions.
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
        """
        Checks if the formula string is valid for parsing:
        if it's not contains any characters that are not
        allowed.
        """
        search = compile(self.allowed_symbols).findall
        if search(self.formula):
            return search(self.formula)
        else:
            return True

    def are_brackets_balanced(self) -> bool:
        """
        Check if all sort of brackets come in pairs.
        """
        c = Counter(self.formula)
        bracket_counter = c["["] == c["]"] and c["{"] == c["}"] and c["("] == c[")"]

        return bracket_counter

    def is_adduct_one(self) -> bool:
        """
        Check if there is only one adduct in formula
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

    def are_atoms_legal(self, parsed) -> None:
        """
        Checks if all parsed atom are belong to
        periodic table.
        """
        for atom in list(parsed.keys()):
            if atom not in self.list_of_atoms:
                raise NoSuchAtom("No atom %s in the periodic table!" % atom)
        return

    def parse_formula(self) -> dict:
        """
        Parse the formula and return a dict with occurences of each atom.
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
