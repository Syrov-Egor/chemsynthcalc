import re
from .formula import ChemicalFormula


class ChemicalFormulaParser(ChemicalFormula):

    def _dictify(self, tuples: list[tuple[str, ...]]) -> dict[str, float]:
        """Transform list of tuples to a dict of atoms.

        Arguments:
            tuples (list): list of tuples of atoms

        Returns:
            dict: dictionary of atoms
        """
        result: dict[str, float] = dict()
        for atom, n, _, _ in tuples:
            try:
                result[atom] += float(n or 1)
            except KeyError:
                result[atom] = float(n or 1)
        return result

    def _fuse(
        self, mol1: dict[str, float], mol2: dict[str, float], weight: float = 1.0
    ) -> dict[str, float]:
        """Fuse 2 dicts representing molecules.

        Arguments:
            mol1 (dict): dict of atoms 1
            mol2 (dict): dict of atoms 2
            weight (int): weight

        Returns:
            dict: a new fused dict
        """

        fused_set: set[str] = set(mol1) | set(mol2)
        fused_dict: dict[str, float] = {
            atom: (mol1.get(atom, 0) + mol2.get(atom, 0)) * weight for atom in fused_set
        }

        return fused_dict

    def _parse(self, formula: str) -> tuple[dict[str, float], int]:
        # TODO if token in self.adduct_symbol

        token_list: list[str] = []
        mol: dict[str, float] = {}
        i: int = 0

        while i < len(formula):
            token: str = formula[i]

            if token in self.closer_brackets:
                coefficient_match: re.Match[str] | None = re.match(
                    self.coefficient_regex, formula[i + 1 :]
                )
                if coefficient_match != None:
                    weight: float = float(coefficient_match.group(0))
                    i += len(coefficient_match.group(0))
                else:
                    weight = 1.0
                submol: dict[str, float] = self._dictify(
                    re.findall(self.atom_and_coefficient_regex, "".join(token_list))
                )
                return self._fuse(mol, submol, weight), i

            elif token in self.opener_brackets:
                recursive_dive: tuple[dict[str, float], int] = self._parse(
                    formula[i + 1 :]
                )
                submol = recursive_dive[0]
                lenght: int = recursive_dive[1]
                mol = self._fuse(mol, submol)
                i += lenght + 1

            else:
                token_list.append(token)

            i += 1

        extract_from_tokens: list[tuple[str, ...]] = re.findall(
            self.atom_and_coefficient_regex, "".join(token_list)
        )
        fused_dict: dict[str, float] = self._fuse(
            mol, self._dictify(extract_from_tokens)
        )
        return fused_dict, i

    def parse_formula(self) -> dict[str, float]:
        return self._parse(self.formula)[0]
