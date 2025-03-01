import re

from .formula import Formula


class ChemicalFormulaParser(Formula):
    """
    Parser of chemical formulas.

    Methods of this class take string of compound's chemical formula
    and turn it into a dict of atoms as keys and their coefficients as values.
    """

    def _dictify(self, tuples: list[tuple[str, ...]]) -> dict[str, float]:
        """
        Transform list of tuples to a dict of atoms.

        Parameters:
            tuples (list[tuple[str, ...]]): List of tuples of atoms

        Returns:
            Dictionary of atoms and they quantities
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
        """Fuse together 2 dicts representing molecules.

        Parameters:
            mol1 (dict[str, float]): Dict of atoms 1
            mol2 (dict[str, float]): Dict of atoms 2
            weight (float): Weight

        Returns:
            A new fused dict
        """

        fused_set: set[str] = set(mol1) | set(mol2)
        fused_dict: dict[str, float] = {
            atom: (mol1.get(atom, 0) + mol2.get(atom, 0)) * weight for atom in fused_set
        }

        return fused_dict

    def _parse(self, formula: str) -> tuple[dict[str, float], int]:
        """
        Parse the formula string

        Recurse on opening brackets to parse the subpart and
        return on closing ones because it is the end of said subpart.
        Formula is the argument of this method due to the complications
        of self. Constructions in recursive functions.

        Parameters:
            formula (str): Formula string

        Returns:
            A tuple of the molecule dict and length of parsed part
        """
        token_list: list[str] = []
        mol: dict[str, float] = {}
        i: int = 0

        while i < len(formula):
            token: str = formula[i]

            if token in self.adduct_symbols:
                coefficient_match: re.Match[str] | None = re.match(
                    self.coefficient_regex, formula[i + 1 :]
                )
                if coefficient_match and coefficient_match.group(0) != "":
                    weight: float = float(coefficient_match.group(0))
                    i += len(coefficient_match.group(0))
                else:
                    weight = 1.0
                recursive_dive: tuple[dict[str, float], int] = self._parse(
                    f"({formula[i + 1 :]}){weight}"
                )
                submol = recursive_dive[0]
                lenght: int = recursive_dive[1]
                mol = self._fuse(mol, submol)
                i += lenght + 1

            elif token in self.closer_brackets:
                coefficient_match: re.Match[str] | None = re.match(
                    self.coefficient_regex, formula[i + 1 :]
                )
                if coefficient_match and coefficient_match.group(0) != "":
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

    def _order_output_dict(self, parsed: dict[str, float]) -> dict[str, float]:
        """
        Arranges the unparsed formula in the order in which the chemical
        elements appear in it.

        Parameters:
            parsed (dict[str, float]): A formula parsed by _parse

        Returns:
            An ordered dictionary
        """
        atoms_list: list[str] = re.findall(self.atom_regex, self.formula)
        weights: list[float] = []
        for atom in atoms_list:
            weights.append(parsed[atom])
        return dict(zip(atoms_list, weights))

    def parse_formula(self) -> dict[str, float]:
        """
        Parsing and ordering of formula

        Returns:
            Parsed formula
        """
        parsed = self._parse(self.formula)[0]
        return self._order_output_dict(parsed)
