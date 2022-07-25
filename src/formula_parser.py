from re import findall, match
from collections import Counter
from .periodic_table import  periodic_table

class ChemicalFormulaParser():
    '''
    Parser of chemical formulas. Methods of this class take string of compound chemical formula
    and turn it to dict of atoms as keys and their coefficients as values. For example:
    "H2O" -> {'O': 1, 'H': 2}
    "C2H5OH" -> {'H': 6.0, 'O': 1.0, 'C': 2.0}
    "(K0.6Na0.4)2SO4(H2O)7" -> {'O': 11.0, 'K': 1.2, 'S': 1.0, 'Na': 0.8, 'H': 14.0}
    Atoms in the output dict are in the random order due to the nature of dicts in Python
    '''
    def __init__(self, formula:str) -> None:
        self.atom_regex = '([A-Z][a-z]*)((\d+(\.\d+)?)*)'
        self.opener_brackets = '({['
        self.closer_brackets = ')}]'
        self.formula = formula
        self.list_of_atoms = [x[0] for x in periodic_table]
    
    def __dictify(self, tuples:tuple) -> dict:
        '''
        Transform tuples of tuples to a dict of atoms.
        '''
        res = dict()
        for atom, n, m, k in tuples:
            try:
                res[atom] += float(n or 1)
            except KeyError:
                res[atom] = float(n or 1)
        return res
        
    def __fuse(self, mol1:dict, mol2:dict, w:int=1) -> dict:
        '''
        Fuse 2 dicts representing molecules. Return a new dict.
        '''
        return {atom: (mol1.get(atom, 0) + mol2.get(atom, 0)) * w for atom in set(mol1) | set(mol2)}

    def __parse(self, formula:str) -> dict:
        '''
        Return the molecule dict and length of parsed part.
        Recurse on opening brackets to parse the subpart and
        return on closing ones because it is the end of said subpart.
        Formula is the argument of this method due to the complications 
        of self. constructions in recursive functions.
        '''
        q = []
        mol = {}
        i = 0

        while i < len(formula):
            # Using a classic loop allow for manipulating the cursor
            token = formula[i]

            if token in self.closer_brackets:
                # Check for an index for this part
                m = match('(\d+(\.\d+)?)*', formula[i+1:])
                if m:
                    weight = float(m.group(0))
                    i += len(m.group(0))
                else:
                    weight = 1
                submol = self.__dictify(findall(self.atom_regex, ''.join(q)))
                return self.__fuse(mol, submol, weight), i

            elif token in self.opener_brackets:
                submol, l = self.__parse(formula[i+1:])
                mol = self.__fuse(mol, submol)
                # skip the already read submol
                i += l + 1
            else:
                q.append(token)

            i+=1

        # Fuse in all that's left at base level
        return self.__fuse(mol, self.__dictify(findall(self.atom_regex, ''.join(q)))), i

    def is_brackets_balanced(self) -> bool:
        '''
        Check if all sort of brackets come in pairs
        '''
        c = Counter(self.formula)
        return c['['] == c[']'] and c['{'] == c['}'] and c['('] == c[')']
    
    def are_atoms_legal(self, parsed) -> bool:
        for atom in list(parsed.keys()):
            if atom not in self.list_of_atoms:
                raise ValueError("No atom %s in the periodic table!" % atom)
        return
            

    def parse_formula(self) -> dict:
        '''
        Parse the formula and return a dict with occurences of each atom.
        '''
        if not self.is_brackets_balanced():
            raise ValueError("The brackers are not balanced ![{]$[&?)]}!]")
        parsed = self.__parse(self.formula)[0]
        self.are_atoms_legal(parsed)
        return parsed