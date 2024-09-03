import re
from collections import Counter

from .periodic_table import PeriodicTable
from .formula import Formula
from .chem_errors import (
    EmptyFormula,
    InvalidCharacter,
    NoSuchAtom,
    MoreThanOneAdduct,
    BracketsNotPaired,
)


class FormulaValidator(Formula):

    def _check_empty_formula(self) -> bool:
        return self.formula == ""

    def _invalid_charachers(self) -> list[str]:
        return re.compile(self.allowed_symbols).findall(self.formula)

    def _invalid_atoms(self) -> list[str]:
        atoms_list: list[str] = re.findall(self.atom_regex, self.formula)
        invalid: list[str] = []
        new_string: str = self.formula
        #a trick to get Cl first of C in cases like CCl4
        atoms_list = sorted(list(set(atoms_list)), key=len, reverse=True)
        for atom in atoms_list:
            if atom not in PeriodicTable().atoms:
                invalid.append(atom)
            new_string = new_string.replace(atom, "")
        found_leftovers: list[str] = re.findall(r"[a-z]", new_string)
        invalid.extend(found_leftovers)
        return invalid

    def _bracket_balance(self) -> bool:
        c: Counter[str] = Counter(self.formula)
        for i in range(len(self.opener_brackets)):
            if c[self.opener_brackets[i]] != c[self.closer_brackets[i]]:
                return False
        return True

    def _num_of_adducts(self) -> int:
        c: Counter[str] = Counter(self.formula)
        i: int = 0
        for adduct in self.adduct_symbols:
            i += c[adduct]
        return i

    def validate_formula(self) -> bool:
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
