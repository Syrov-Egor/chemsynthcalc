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
        return [atom for atom in atoms_list if atom not in PeriodicTable().atoms]

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
                f"No atom(s) {self._invalid_atoms()} in the periodic table!"
            )
        elif not self._bracket_balance():
            raise BracketsNotPaired(
                f"The brackets {self.opener_brackets} {self.closer_brackets} are not balanced!"
            )
        elif self._num_of_adducts() > 1:
            raise MoreThanOneAdduct(
                f"More than one adduct {self.adduct_symbols} in the formula"
            )
        return True
