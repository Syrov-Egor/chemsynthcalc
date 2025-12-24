import re
from collections import Counter

from .chem_errors import (
    BracketsNotPaired,
    EmptyFormula,
    InvalidCharacter,
    MoreThanOneAdduct,
    NoSuchAtom,
)
from .formula import Formula
from .periodic_table import ATOMS


class FormulaValidator(Formula):
    """
    Methods of this class validate the initial input formula.
    """

    def _check_empty_formula(self) -> bool:
        """
        Checks if formula is an empty string.
        """
        return self.formula == ""

    def _invalid_charachers(self) -> list[str]:
        """
        Checks if formula contains invalid characters.

        Returns:
            List of invalid characters
        """
        return re.compile(self.allowed_symbols).findall(self.formula)

    def _invalid_atoms(self) -> list[str]:
        """
        Checks whether the formula contains atoms that
        are not in the periodic system.

        Returns:
            List of invalid atoms
        """
        atoms_list: list[str] = re.findall(self.atom_regex, self.formula)
        invalid: list[str] = []
        new_string: str = self.formula
        # a trick to get Cl first of C in cases like CCl4
        atoms_list = sorted(list(set(atoms_list)), key=len, reverse=True)
        for atom in atoms_list:
            if atom not in ATOMS:
                invalid.append(atom)
            new_string = new_string.replace(atom, "")
        found_leftovers: list[str] = re.findall(r"[a-z]", new_string)
        invalid.extend(found_leftovers)
        return invalid

    def _bracket_balance(self) -> bool:
        """
        Checks whether all of the brackets come in pairs.
        """
        c: Counter[str] = Counter(self.formula)
        for i in range(len(self.opener_brackets)):
            if c[self.opener_brackets[i]] != c[self.closer_brackets[i]]:
                return False
        return True

    def _num_of_adducts(self) -> int:
        """
        Counts a number of adduct symbols
        (listed in [Formula base class][chemsynthcalc.formula.Formula]).

        Returns:
            A number of adduct symbols
        """
        c: Counter[str] = Counter(self.formula)
        i: int = 0
        for adduct in self.adduct_symbols:
            i += c[adduct]
        return i

    def validate_formula(self) -> bool:
        """
        Validation of the formula string.
        Calls the private methods of this class in order.

        Raise:
            [EmptyFormula][chemsynthcalc.chem_errors.EmptyFormula] if the formula is an empty string. <br />
            [InvalidCharacter][chemsynthcalc.chem_errors.InvalidCharacter] if there is an invalid character(s) in the string. <br />
            [NoSuchAtom][chemsynthcalc.chem_errors.NoSuchAtom] if there is an invalid atom(s) in the string. <br />
            [BracketsNotPaired][chemsynthcalc.chem_errors.BracketsNotPaired] if the brackets are not in pairs. <br />
            [MoreThanOneAdduct][chemsynthcalc.chem_errors.MoreThanOneAdduct] if there are more than 1 adduct symbols in the string. <br />

        Returns:
            True if all the checks are OK
        """
        if self._check_empty_formula():
            raise EmptyFormula
        elif self._invalid_charachers():
            raise InvalidCharacter(
                f"Invalid character(s) {self._invalid_charachers()} in the formula {self.formula}"
            )
        elif self._invalid_atoms():
            raise NoSuchAtom(
                f"The formula {self.formula} contains atom {self._invalid_atoms()} which is not in the periodic table"
            )
        elif not self._bracket_balance():
            raise BracketsNotPaired(
                f"The brackets {self.opener_brackets} {self.closer_brackets} are not balanced the formula {self.formula}!"
            )
        elif self._num_of_adducts() > 1:
            raise MoreThanOneAdduct(
                f"More than one adduct {self.adduct_symbols} in the formula {self.formula}"
            )
        return True
