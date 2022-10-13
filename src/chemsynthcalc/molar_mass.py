from .periodic_table import periodic_table
from .formula_parser import ChemicalFormulaParser


class MolarMassCalculation:
    """Class for the calculation of molar masses and percentages of compounds.
    
    Compounds should be parsed by :class:`chemsynthcalc.formula_parser.ChemicalFormulaParser` first.

    Arguments:
        parsed_formula (dict): formula parsed by :class:`chemsynthcalc.formula_parser.ChemicalFormulaParser`

    Attributes:
        atomicWeight (list): list of (atom, weight) tuple formed from :mod:`chemsynthcalc.periodic_table`
        list_of_atoms (list): list of atoms formed from :mod:`chemsynthcalc.periodic_table`
    """

    def __init__(self, parsed_formula: dict) -> None:
        self.atomicWeight = [(x[0], x[1]) for x in periodic_table]
        self.parsed_formula = parsed_formula
        self.list_of_atoms = [x[0] for x in periodic_table]

    def calculate_molar_mass(self) -> float:
        """Calculation of the molar mass of compound from the atomic masses of atoms in a parsed formula.
        
        Atomic masses are taken from :mod:`chemsynthcalc.periodic_table` file.

        Returns:
            float: molar mass (in g/mol)

        Examples:
            >>> MolarMassCalculation({'H':2, 'O':1}).calculate_molar_mass()
            18.015
            >>> MolarMassCalculation({'C':2, 'H':6, 'O':1}).calculate_molar_mass()
            46.069
        """
        atomic_masses_list = []
        for atom, n in self.parsed_formula.items():
            atom_index = self.list_of_atoms.index(atom)
            atomic_masses_list.append(self.atomicWeight[atom_index][1] * n)
        return sum(atomic_masses_list)

    def calculate_mass_percent(self) -> dict:
        """Calculation of mass percents of atoms in parsed formula.

        Returns:
            dict: mass percentages of atoms in formula
        
        Examples:
            >>> MolarMassCalculation({'H':2, 'O':1}).calculate_mass_percent()
            {'H': 11.19067443796836, 'O': 88.80932556203163}
            >>> MolarMassCalculation({'C':2, 'H':6, 'O':1}).calculate_mass_percent()
            {'C': 52.14352384466777, 'H': 13.12813388612733, 'O': 34.72834226920489}
        """
        atomic_masses_list = []
        for atom, n in self.parsed_formula.items():
            atom_index = self.list_of_atoms.index(atom)
            atomic_masses_list.append(self.atomicWeight[atom_index][1] * n)
        molar_mass = self.calculate_molar_mass()
        percents = [atomic / molar_mass * 100 for atomic in atomic_masses_list]
        return dict(zip(self.parsed_formula.keys(), percents))

    def calculate_atomic_percent(self) -> dict:
        """Calculation of atomic percents of atoms in the parsed formula.

        Returns:
            dict: atomic percentages of atoms in formula
        
        Examples:
            >>> MolarMassCalculation({'H':2, 'O':1}).calculate_atomic_percent()
            {'H': 66.66666666666666, 'O': 33.33333333333333
            >>> MolarMassCalculation({'C':2, 'H':6, 'O':1}).calculate_atomic_percent()
            {'C': 22.22222222222222, 'H': 66.66666666666666, 'O': 11.11111111111111}
        """
        values = list(self.parsed_formula.values())
        total = sum(values)
        atomic = [i / total * 100 for i in values]
        return dict(zip(self.parsed_formula.keys(), atomic))

    def calculate_oxide_percent(self) -> dict:
        """Calculation of oxide percents in parsed formula.

        Calculation of oxide percents in parsed formula from the types of oxide 
        declared in the periodic table file. This type of data
        is mostly used in XRF spectrometry and mineralogy. The oxide
        percents are calculated by finding the convertion factor between element
        and its respective oxide (https://www.geol.umd.edu/~piccoli/probe/molweight.html)
        and normalizing the total sum to 100%. One can change the oxide type
        for certain elements in the :mod:`chemsynthcalc.periodic_table` file. Theoretically, this function
        should work for other types of binary compound (sulfides, fluorides etc.)
        or even salts, however, modification of this function is required
        (for instance, in case of binary compound, removing X atom
        from the list of future compounds should have X as an argument of this function)

        Returns:
            dict: oxide percentages of oxides in formula
        
        Examples:
            >>> MolarMassCalculation({'C':2, 'H':6, 'O':1}).calculate_oxide_percent()
            {'CO2': 61.9570190690046, 'H2O': 38.04298093099541}
            >>> MolarMassCalculation({'Ba':1, 'Ti':1, 'O':3}).calculate_oxide_percent()
            {'BaO': 65.7516917244869, 'TiO2': 34.24830827551309}
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
