import sys
import time
import json
import numpy as np


class FormulaOutput:
    """Outputs for :class:`chemsynthcalc.chemical_formula.ChemicalFormula`

    A collection of methods for outputing results of
    :class:`chemsynthcalc.chemical_formula.ChemicalFormula`. Default filenames of
    output files are "CSC_formula" with formula
    and nanosecond since the epoch.
    
    Arguments:
        chemical_formula_output (dict): output of :class:`chemsynthcalc.chemical_formula.ChemicalFormula`
    """

    def __init__(self, chemical_formula_output: dict) -> None:
        self.output: dict = chemical_formula_output
        self.original_stdout = sys.stdout

    def __print_stream(self, print_rounding_order: int) -> None:
        """Prints of results that can be streamed into stdout.

        Arguments:
            print_rounding_order (int): print precision
        
        Returns:
            None
        """
        print("formula:", self.output.get("formula"))
        print("parsed formula:", self.output.get("parsed formula"))
        print("molar mass:", round(self.output.get("molar mass"), print_rounding_order))
        print(
            "mass percent:",
            {
                k: round(v, print_rounding_order)
                for k, v in self.output.get("mass percent").items()
            },
        )
        print(
            "atomic percent:",
            {
                k: round(v, print_rounding_order)
                for k, v in self.output.get("atomic percent").items()
            },
        )
        print(
            "oxide percent:",
            {
                k: round(v, print_rounding_order)
                for k, v in self.output.get("oxide percent").items()
            },
        )
        return

    def print_results(self, print_rounding_order: int) -> None:
        """Print results in the terminal.
        
        Arguments:
            print_rounding_order (int): print precision
                
        Returns:
            None
        """
        sys.stdout = self.original_stdout
        self.__print_stream(print_rounding_order)
        return

    def export_to_txt(self, filename: str, print_rounding_order: int) -> None:
        """Prints results into the txt file and saves them.
        
        Arguments:
            filename (str): filename string (should end with .txt)
            print_rounding_order (int): print precision
                
        Returns:
            None
        """
        if filename == "default":
            filename = "CSC_formula_%s_%s.txt" % (
                self.output.get("formula:"),
                time.time_ns(),
            )

        with open(filename, "w") as file:
            sys.stdout = file
            self.__print_stream(print_rounding_order)
        sys.stdout = self.original_stdout
        return

    def dump_to_json(self, print_rounding_order: int) -> str:
        """JSON serialization.
        
        Arguments:
            print_rounding_order (int): print precision

        Returns:
            str: A JSON-type object of results output.
        """
        mod_output = self.output.copy()
        mod_output.update(
            {
                "mass percent": {
                    k: round(v, print_rounding_order)
                    for k, v in self.output.get("mass percent").items()
                }
            }
        )
        mod_output.update(
            {
                "atomic percent": {
                    k: round(v, print_rounding_order)
                    for k, v in self.output.get("atomic percent").items()
                }
            }
        )
        mod_output.update(
            {
                "oxide percent": {
                    k: round(v, print_rounding_order)
                    for k, v in self.output.get("oxide percent").items()
                }
            }
        )
        return json.dumps(mod_output, ensure_ascii=False)

    def export_to_json(self, filename: str, print_rounding_order: int) -> None:
        """Dump the output dict into an JSON file.
 
        Arguments:
            filename (str): filename string (should end with .json)
            print_rounding_order (int): print precision
                
        Returns:
            None
        """
        if filename == "default":
            filename = "CSC_formula_%s_%s.json" % (
                self.output.get("formula"),
                time.time_ns(),
            )

        with open(filename, "w", encoding="utf-8") as file:
            json_obj = json.loads(self.dump_to_json(print_rounding_order))
            json.dump(json_obj, file, ensure_ascii=False)
        return


class ReactionOutput:
    """Outputs for :class:`chemsynthcalc.chemical_reaction.ChemicalReaction`
    
    A collection of methods to output the results of
    :class:`chemsynthcalc.chemical_reaction.ChemicalReaction`.The default 
    filenames of the output files are "CSC_reaction" with target compound and nanoseconds since the epoch.

    Arguments:
        chemical_formula_output (dict): output of :class:`chemsynthcalc.chemical_reaction.ChemicalReaction`
    """

    def __init__(self, chemical_reaction_output: dict) -> None:
        self.output: dict = chemical_reaction_output
        self.original_stdout = sys.stdout

    def __print_stream(self, print_rounding_order: int) -> None:
        """Prints of results that can be streamed into stdout.

        Arguments:
            print_rounding_order (int): print precision
        
        Returns:
            None
        """
        print("initial reaction:", self.output.get("initial reaction"))
        print("reaction matrix:")
        print(self.output.get("reaction matrix"))
        print("mode:", self.output.get("mode"))
        print("coefficients:", self.output.get("coefficients"))
        print("normalized coefficients:", self.output.get("normalized coefficients"))
        print("balanced by algorithm: %s " % self.output.get("algorithm"))
        print("is balanced:", self.output.get("is balanced"))
        print("final reaction:", self.output.get("final reaction"))
        print(
            "final reaction normalized:", self.output.get("final reaction normalized")
        )
        print("target:", self.output.get("target"))
        formulas = self.output.get("formulas")
        molar_masses = self.output.get("molar masses")
        masses = self.output.get("masses")
        for i, formula in enumerate(formulas):
            print(
                "%s: M = %s g/mol, m = %s g"
                % (
                    formula,
                    "%.{0}f".format(print_rounding_order)
                    % round(molar_masses[i], print_rounding_order),
                    "%.{0}f".format(print_rounding_order)
                    % round(masses[i], print_rounding_order),
                )
            )
        return

    def print_results(self, print_rounding_order: int) -> None:
        """Print results in the terminal.
        
        Arguments:
            print_rounding_order (int): print precision
                
        Returns:
            None
        """
        sys.stdout = self.original_stdout
        self.__print_stream(print_rounding_order)
        return

    def export_to_txt(self, filename: str, print_rounding_order: int) -> None:
        """Prints results into the txt file and saves it.
        
        Arguments:
            filename (str): filename string (should end with .txt)
            print_rounding_order (int): print precision
                
        Returns:
            None
        """
        if filename == "default":
            filename = "CSC_reaction_%s_%s.txt" % (
                self.output.get("target"),
                time.time_ns(),
            )
        with open(filename, "w") as file:
            sys.stdout = file
            self.__print_stream(print_rounding_order)
        sys.stdout = self.original_stdout
        return

    def dump_to_json(self, print_rounding_order: int) -> str:
        """JSON serialization.
        
        Arguments:
            print_rounding_order (int): print precision

        Returns:
            str: A JSON-type object of results output.
        """
        mod_output = self.output.copy()
        mod_output.update(
            {"reaction matrix": np.array2string(self.output.get("reaction matrix"))}
        )
        str_formulas = [str(formula) for formula in self.output.get("formulas")]
        mod_output.update({"formulas": str_formulas})
        mod_output.update({"target": str(self.output.get("target"))})
        mod_output.update(
            {
                "molar masses": [
                    round(mass, print_rounding_order)
                    for mass in self.output.get("molar masses")
                ]
            }
        )
        mod_output.update(
            {
                "masses": [
                    round(mass, print_rounding_order)
                    for mass in self.output.get("masses")
                ]
            }
        )
        return json.dumps(mod_output, ensure_ascii=False)

    def export_to_json(self, filename: str, print_rounding_order: int) -> None:
        """Dump output dict into JSON flie.
 
        Arguments:
            filename (str): filename string (should end with .json)
            print_rounding_order (int): print precision
                
        Returns:
            None
        """
        if filename == "default":
            filename = "CSC_reaction_%s_%s.json" % (
                self.output.get("target"),
                time.time_ns(),
            )

        with open(filename, "w", encoding="utf-8") as file:
            json_obj = json.loads(self.dump_to_json(print_rounding_order))
            json.dump(json_obj, file, ensure_ascii=False)
        return