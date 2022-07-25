from .periodic_table import periodic_table
import numpy as np

class ChemicalReactionMatrix():
    '''
    Class that creates a reaction matrix from 
    parced formulas of the reaction string
    '''
    def __init__(self, parsed_all:list) -> None:
        self.parsed_all = parsed_all
        self.merged_dict = {k: v for d in self.parsed_all for k, v in d.items()}
        self.electronegativity = [(x[0], x[2]) for x in periodic_table]
        self.electronegativity_pattern = [atom[0] for atom in
            sorted(self.electronegativity, key=lambda tup: tup[1])]
        self.list_of_atoms = [x[0] for x in periodic_table]

    def create_reaction_matrix(self) -> np.array:
        elements = sorted(list(self.merged_dict.keys()),
        key=lambda x: self.electronegativity_pattern.index(x))
        elements_by_number = []
        for atom in elements:
            elements_by_number.append(self.list_of_atoms.index(atom)+1)
        reaction_matrix = [[elements_by_number[i]] for i, atom in enumerate(elements)]
        for array_index, array in enumerate(reaction_matrix):
            for parsed_reactant in self.parsed_all:
                if elements[array_index] in parsed_reactant.keys():
                    array.append(parsed_reactant.get(elements[array_index]))
                else:
                    array.append(0)
        return np.array([array[1:] for array in reaction_matrix])