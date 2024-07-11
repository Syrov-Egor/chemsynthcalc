import re
from .periodic_table import PeriodicTable

class ChemicalFormulaParser:
    def __init__(self, formula: str) -> None:
        self.atom_regex: str = r"([A-Z][a-z]*)"
        self.coefficient_regex: str = r"((\d+(\.\d+)?)*)"
        self.atom_and_coefficient_regex: str = self.atom_regex + self.coefficient_regex
        self.opener_brackets: str = "({["
        self.closer_brackets: str = ")}]"
        self.adduct_symbols: str = "*·•"
        self.allowed_symbols: str = r"[^A-Za-z0-9.({[)}\]*·•]"

        self.formula: str = formula
        self.atoms: set[str] = PeriodicTable().atoms

    def _dictify(self, tuples: tuple[tuple]) -> dict:
        """Transform tuples of tuples to a dict of atoms.

        Arguments:
            tuples (tuple): tuple of tuples of atoms

        Returns:
            dict: dictionary of atoms
        """
        result: dict = dict()
        for atom, n, m, k in tuples:
            try:
                result[atom] += float(n or 1)
            except KeyError:
                result[atom] = float(n or 1)
        return result

    def _fuse(self, mol1: dict, mol2: dict, w: int = 1) -> dict:
        """Fuse 2 dicts representing molecules. 
        
        Arguments:
            mol1 (dict): dict of atoms 1
            mol2 (dict): dict of atoms 2
            w (int): weight

        Returns:
            dict: a new fuse dict
        """

        fused_set: set = set(mol1) | set(mol2)
        fused_dict: dict = {atom: (mol1.get(atom, 0) + mol2.get(atom, 0)) * w for atom in fused_set}

        return fused_dict

    def _parse(self, formula: str) -> dict:
        #TODO elif token in self.adduct_symbol

        token_list: list = []
        mol: dict = {}
        i: int = 0

        while i < len(formula):
            token: str = formula[i]

            if token in self.closer_brackets:
                coefficient_match: str = re.match(self.coefficient_regex, formula[i + 1 :]).group(0)
                if coefficient_match != "":
                    weight: float = float(coefficient_match)
                    i += len(coefficient_match)
                else:
                    weight = 1.0
                submol: dict = self._dictify(re.findall(self.atom_and_coefficient_regex, "".join(token_list)))

            elif token in self.opener_brackets:
                recursive_dive: tuple[dict, int] = self._parse(formula[i + 1 :])
                submol: dict = recursive_dive[0]
                lenght: int = recursive_dive[1]
                mol = self._fuse(mol, submol)
                i += lenght + 1
            
            else:
                token_list.append(token)

            i += 1

        return token_list

    def parse_formula(self) -> dict:
        output: dict = {}
        output = self._parse(self.formula)
        return output