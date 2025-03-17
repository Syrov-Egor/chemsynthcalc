from typing import NamedTuple

from .periodic_table import PeriodicTable
from .formula_parser import ChemicalFormulaParser


class Oxide(NamedTuple):
    """
    A named tuple to represent oxide properties: a first atom
    (usually metal), full oxide compound string and mass percent
    of the first atom in the initial formula.
    """

    atom: str
    label: str
    mass_percent: float


class MolarMassCalculation:
    """
    Class for the calculation of molar masses and percentages of a compound.

    Compound should be parsed by
    [ChemicalFormulaParser][chemsynthcalc.formula_parser.ChemicalFormulaParser] first.

    Parameters:
        parsed_formula (dict): Formula parsed by [ChemicalFormulaParser][chemsynthcalc.formula_parser.ChemicalFormulaParser]

    Attributes:
        p_table (dict[str, Atom]): Periodic table of elements
    """

    def __init__(self, parsed_formula: dict[str, float]) -> None:
        self.parsed_formula: dict[str, float] = parsed_formula
        self.p_table = PeriodicTable().p_table

    def _calculate_atomic_masses(self) -> list[float]:
        """
        Calculation of the molar masses of
        all atoms in a parsed formula.

        Returns:
            List of atomic masses multiplied by the number of corresponding atoms
        """
        masses: list[float] = []
        for atom, weight in self.parsed_formula.items():
            atom_mass: float = self.p_table[atom].atomic_weight
            masses.append(atom_mass * weight)
        return masses

    def calculate_molar_mass(self) -> float:
        """
        Calculation of the molar mass of compound from
        the atomic masses of atoms in a parsed formula.

        Returns:
            Molar mass (in g/mol)

        Examples:
            >>> MolarMassCalculation({'H':2, 'O':1}).calculate_molar_mass()
            18.015
            >>> MolarMassCalculation({'C':2, 'H':6, 'O':1}).calculate_molar_mass()
            46.069
        """
        return sum(self._calculate_atomic_masses())

    def calculate_mass_percent(self) -> dict[str, float]:
        """
        Calculation of mass percents of atoms in parsed formula.

        Returns:
            Mass percentages of atoms in the formula

        Examples:
            >>> MolarMassCalculation({'H':2, 'O':1}).calculate_mass_percent()
            {'H': 11.19067443796836, 'O': 88.80932556203163}
            >>> MolarMassCalculation({'C':2, 'H':6, 'O':1}).calculate_mass_percent()
            {'C': 52.14352384466777, 'H': 13.12813388612733, 'O': 34.72834226920489}
        """
        atomic_masses: list[float] = self._calculate_atomic_masses()
        molar_mass: float = self.calculate_molar_mass()
        percents: list[float] = [atomic / molar_mass * 100 for atomic in atomic_masses]
        return dict(zip(self.parsed_formula.keys(), percents))

    def calculate_atomic_percent(self) -> dict[str, float]:
        """
        Calculation of atomic percents of atoms in the parsed formula.

        Returns:
            Atomic percentages of atoms in the formula

        Examples:
            >>> MolarMassCalculation({'H':2, 'O':1}).calculate_atomic_percent()
            {'H': 66.66666666666666, 'O': 33.33333333333333
            >>> MolarMassCalculation({'C':2, 'H':6, 'O':1}).calculate_atomic_percent()
            {'C': 22.22222222222222, 'H': 66.66666666666666, 'O': 11.11111111111111}
        """
        values: list[float] = list(self.parsed_formula.values())
        atomic: list[float] = [value / sum(values) * 100 for value in values]
        return dict(zip(self.parsed_formula.keys(), atomic))

    def _custom_oxides_input(self, *args: str) -> list[Oxide]:
        """
        Checks if passed non-default oxide formulas can be applied to any
        atom in parsed formula. If so, replaces it with said non-default oxide.
        If not, the default oxide is chosen.

        Parameters:
            *args (tuple[str, ...]): An arbitrary number of non-default oxide formulas

        Returns:
            A list of Oxide objects

        Raise:
            ValueError if compound is not binary or second element is not oxygen
        """
        first_atoms: list[str] = []
        for c_oxide in args:
            parsed_oxide = list(ChemicalFormulaParser(c_oxide).parse_formula().keys())

            if len(parsed_oxide) > 2:
                raise ValueError("Only binary compounds can be considered as input")

            elif parsed_oxide[1] != "O":
                raise ValueError("Only oxides can be considered as input")

            first_atoms.append(parsed_oxide[0])

        custom_oxides = dict(zip(first_atoms, args))
        mass_percents: list[float] = list(self.calculate_mass_percent().values())
        oxides: list[Oxide] = []
        for i, atom in enumerate(self.parsed_formula.keys()):
            if atom != "O":
                if atom in custom_oxides.keys():
                    label = custom_oxides[atom]
                else:
                    label = self.p_table[atom].default_oxide
                oxides.append(Oxide(atom, label, mass_percents[i]))

        return oxides

    def calculate_oxide_percent(self, *args: str) -> dict[str, float]:
        """
        Calculation of oxide percents in parsed formula.

        Calculation of oxide percents in parsed formula from the types of oxide
        declared in the periodic table file. This type of data
        is mostly used in XRF spectrometry and mineralogy. The oxide
        percents are calculated by finding the [convertion factor between element
        and its respective oxide](https://www.geol.umd.edu/~piccoli/probe/molweight.html)
        and normalizing the total sum to 100%. One can change the oxide type
        for certain elements in the [periodic_table.py][chemsynthcalc.periodic_table] file. Theoretically, this function
        should work for other types of binary compound (sulfides, fluorides etc.)
        or even salts, however, modification of this function is required
        (for instance, in case of binary compound, removing X atom
        from the list of future compounds should have X as an argument of this function).

        Parameters:
            *args (tuple[str, ...]): An arbitrary number of non-default oxide formulas

        Returns:
            Percentages of oxides in the formula

        Examples:
            >>> MolarMassCalculation({'C':2, 'H':6, 'O':1}).calculate_oxide_percent()
            {'CO2': 61.9570190690046, 'H2O': 38.04298093099541}
            >>> MolarMassCalculation({'Ba':1, 'Ti':1, 'O':3}).calculate_oxide_percent()
            {'BaO': 65.7516917244869, 'TiO2': 34.24830827551309}
        """
        oxides: list[Oxide] = self._custom_oxides_input(*args)

        oxide_percents: list[float] = []

        for oxide in oxides:
            parsed_oxide: dict[str, float] = ChemicalFormulaParser(
                oxide.label
            ).parse_formula()
            oxide_mass: float = MolarMassCalculation(
                parsed_oxide
            ).calculate_molar_mass()
            atomic_oxide_coef: float = parsed_oxide[oxide.atom]
            atomic_mass: float = self.p_table[oxide.atom].atomic_weight
            conversion_factor: float = oxide_mass / atomic_mass / atomic_oxide_coef
            oxide_percents.append(oxide.mass_percent * conversion_factor)

        normalized_oxide_percents: list[float] = [
            x / sum(oxide_percents) * 100 for x in oxide_percents
        ]
        oxide_labels: list[str] = [oxide.label for oxide in oxides]

        return dict(zip(oxide_labels, normalized_oxide_percents))
