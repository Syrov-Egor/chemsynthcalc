import sys
import time
import json
import numpy as np

class ReactionOutput():
    '''
    A collection of methods for outputing results of
    `ChemicalReaction` class.
    '''
    def __init__(self, chemical_reaction_output:dict) -> None:
        self.output:dict = chemical_reaction_output
        self.original_stdout = sys.stdout
    
    def __print_stream(self, print_rounding_order:int) -> None:
        print("initial reaction:", self.output.get("initial reaction:"))
        print("reaction matrix:")
        print(self.output.get("reaction matrix:"))
        print("mode:", self.output.get("mode:"))
        print("coefficients:", self.output.get("coefficients:"))
        print("normalized coefficients:", self.output.get("normalized coefficients:"))
        print("balanced by algorithm: %s " % self.output.get("algorithm:"))
        print("final reaction:", self.output.get("final reaction:"))
        print("final reaction normalized:", self.output.get("final reaction normalized:"))
        print("target:", self.output.get("target:"))
        formulas = self.output.get("formulas:")
        molar_masses = self.output.get("molar masses:")
        masses = self.output.get("masses:")
        for i, formula in enumerate(formulas):
            print("%s: M = %s g/mol, m = %s g" % (
                formula,
                '%.{0}f'.format(print_rounding_order) % round(molar_masses[i], print_rounding_order),
                '%.{0}f'.format(print_rounding_order) % round(masses[i], print_rounding_order)))
        return
    
    def print_results(self, print_rounding_order:int) -> None:
        '''
        Print results in the terminal.
        '''
        sys.stdout = self.original_stdout
        self.__print_stream(print_rounding_order)
        return
        
    def export_to_txt(self, print_rounding_order:int)  -> None:
        '''
        Print results in the txt file and saves it.
        '''
        filename = "chemsynthcalc_%s_%s.txt" % (self.output.get("target:"), time.time_ns())
        with open(filename, 'w') as file:
            sys.stdout = file
            self.__print_stream(print_rounding_order)
        sys.stdout = self.original_stdout
        return
    
    def export_to_json(self, print_rounding_order:int)  -> None:
        '''
        Dump output dict into JSON flie.
        '''
        filename = "chemsynthcalc_%s_%s.json" % (self.output.get("target:"), time.time_ns())
        mod_output = self.output.copy()
        mod_output.update({"reaction matrix:": np.array2string(self.output.get("reaction matrix:"))})
        str_formulas = [str(formula) for formula in self.output.get("formulas:")]
        mod_output.update({"formulas:": str_formulas})
        mod_output.update({"target:": str(self.output.get("target:"))})
        mod_output.update({"molar masses:": [round(mass, print_rounding_order) for mass in self.output.get("molar masses:")]})
        mod_output.update({"masses:": [round(mass, print_rounding_order) for mass in self.output.get("masses:")]})
        with open(filename, 'w', encoding='utf-8') as file:
            json.dump(mod_output, file, ensure_ascii=False)
        return