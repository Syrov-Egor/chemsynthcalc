from .periodic_table import periodic_table
from .formula_parser import ChemicalFormulaParser


class MolarMassCalculation:
    """
    Class for calculation of molar masses of compounds.
    Compounds should be parsed by FormulaParser first.

    Parameters:
    * `parsed_formula:dict` - formula parsed by `ChemicalFormulaParser`

    For example:
    ```
    >>>MolarMassCalculation({'C':2, 'H':6, 'O':1}).calculate_molar_mass()
    18.015
    >>>MolarMassCalculation({'C':2, 'H':6, 'O':1}).calculate_molar_mass()
    46.069
    ```
    """

    def __init__(self, parsed_formula: dict) -> None:
        self.atomicWeight = [(x[0], x[1]) for x in periodic_table]
        self.parsed_formula = parsed_formula
        self.list_of_atoms = [x[0] for x in periodic_table]

    def calculate_molar_mass(self) -> float:
        """
        Calculation of molar mass of compound from atomic masses of atoms in parsed formula.
        Atomic masses are taken from periodic_table file.
        """
        atomic_masses_list = []
        for atom, n in self.parsed_formula.items():
            atom_index = self.list_of_atoms.index(atom)
            atomic_masses_list.append(self.atomicWeight[atom_index][1] * n)
        return sum(atomic_masses_list)

    def calculate_mass_percent(self) -> dict:
        """
        Calculation of mass percents of atoms in parsed formula.
        """
        atomic_masses_list = []
        for atom, n in self.parsed_formula.items():
            atom_index = self.list_of_atoms.index(atom)
            atomic_masses_list.append(self.atomicWeight[atom_index][1] * n)
        molar_mass = self.calculate_molar_mass()
        percents = [atomic / molar_mass * 100 for atomic in atomic_masses_list]
        return dict(zip(self.parsed_formula.keys(), percents))

    def calculate_atomic_percent(self) -> dict:
        """
        Calculation of atomic percents of atoms in parsed formula.
        """
        values = list(self.parsed_formula.values())
        total = sum(values)
        atomic = [i / total * 100 for i in values]
        return dict(zip(self.parsed_formula.keys(), atomic))

    def calculate_oxide_percent(self) -> dict:
        """
        Calculation of oxide percents in parsed formula from types
        of oxides declared in periodic table file. This type of data
        is mostly used in XRF spectrometry and mineralogy. The oxide
        percents are calculated by finding the convertion factor between elenent
        and its respective oxide (https://www.geol.umd.edu/~piccoli/probe/molweight.html)
        and normalizing the total sum to 100%. One can change the oxide type
        for certain element in periodic_table file. Theoretically, this function
        should work for other types of binary compound (sulfides, fluorides etc.)
        or even salts, however, modification of this function is required
        (for instance, in case of binary compound, removing X atom
        from list of future compounds should have X as an argument of this function)
        """
        oxide_types = [i[3] for i in periodic_table]
        old_atoms = list(self.parsed_formula.keys())
        old_weights = list(
            MolarMassCalculation(self.parsed_formula).calculate_mass_percent().values()
        )
        # Remove O from list of future oxides
        atoms = []
        weights = []
        if "O" in old_atoms:
            temprow = []
            for j, item in enumerate(old_weights):
                if j != old_atoms.index("O"):
                    temprow.append(item)
            weights.append(temprow)
            for k, item in enumerate(old_atoms):
                if k != old_atoms.index("O"):
                    atoms.append(item)
        else:
            atoms = old_atoms
            weights = [old_weights]
        molar_weights = [
            elements[1]
            for atom in atoms
            for elements in periodic_table
            if atom == elements[0]
        ]
        indxs = [
            periodic_table.index(elements)
            for atom in atoms
            for elements in periodic_table
            if atom == elements[0]
        ]
        oxides = [oxide_types[indx] for indx in indxs]
        parsed_oxides = [
            ChemicalFormulaParser(formula).parse_formula() for formula in oxides
        ]
        molar_masses = [
            MolarMassCalculation(oxide).calculate_molar_mass()
            for oxide in parsed_oxides
        ]
        convertion_factors = []
        for i, atom in enumerate(atoms):
            convertion_factor = (
                molar_masses[i] / molar_weights[i] / parsed_oxides[i].get(atom)
            )
            convertion_factors.append(convertion_factor)
        oxide_percents = []
        for row in weights:
            oxide_percent = [row[x] * convertion_factors[x] for x in range(len(row))]
            oxide_percents.append([x / sum(oxide_percent) * 100 for x in oxide_percent])

        return dict(zip(oxides, oxide_percents[0]))
