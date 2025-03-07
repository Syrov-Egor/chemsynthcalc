import sys
import time
import json

import numpy as np

from .utils import round_dict_content


class ChemicalOutput:
    """
    Methods of this class prepare output from
    [ChemicalFormula][chemsynthcalc.chemical_formula.ChemicalFormula]
    and
    [ChemicalReaction][chemsynthcalc.chemical_reaction.ChemicalReaction]
    objects and output it in different ways.

    Arguments:
        output (dict[str, object]): Output dictionary
        print_precision (int): How many decimal places to print out
        obj (str): Type of object ("formula" or "reaction")

    Attributes:
        rounded_values (dict[str, object]): Output dictionary rounded to print_precision
        original_stdout (TextIO | Any): Default stdout

    Raise:
        ValueError if print_precision <= 0 <br / >
        ValueError if obj is not "ChemicalFormula" or "ChemicalReaction"
    """

    def __init__(
        self, output: dict[str, object], print_precision: int, obj: str
    ) -> None:
        if print_precision > 0:
            self.print_precision: int = print_precision
        else:
            raise ValueError("precision <= 0")

        if obj in {"ChemicalFormula", "ChemicalReaction"}:
            self.obj = obj
        else:
            raise ValueError(f"No object of a class: {obj}")

        self.output: dict[str, object] = output
        self.rounded_values: dict[str, object] = self._round_values()
        self.original_stdout = sys.stdout

    def _round_values(self) -> dict[str, object]:
        """
        Round values of output dictionary to the print_precision.
        Rounding is different depending on the type of value.

        Returns:
            Rounded dictionary
        """
        rounded_dict: dict[str, object] = {}
        for name, value in self.output.items():
            if isinstance(value, float):
                rounded_value = round(value, self.print_precision)
            elif isinstance(value, dict):
                rounded_value = round_dict_content(value, self.print_precision)  # type: ignore
            elif name == "masses":
                rounded_value = [round(v, self.print_precision) for v in value]  # type: ignore
            elif name == "reaction matrix":
                rounded_value = np.array2string(value)  # type: ignore
            else:
                rounded_value = value

            rounded_dict.update({name: rounded_value})  # type: ignore

        return rounded_dict

    def _generate_filename(self, file_type: str) -> str:
        """
        Generates a filename for an output file in the form of:
        "CSC_object type_formula or target_nanosec since the Epoch.txt or json"

        Returns:
            String of a filename
        """
        if self.obj == "ChemicalFormula":
            filename: str = (
                f"CSC_{self.obj}_{self.output.get("formula")}_{time.time_ns()}.{file_type}"
            )
        else:
            filename: str = (
                f"CSC_{self.obj}_{self.output.get("target")}_{time.time_ns()}.{file_type}"
            )

        return filename

    def _print_additional_reaction_results(self) -> None:
        """
        Print output masses in a user-friendly human-readable format.
        """
        for i, formula in enumerate(self.output["formulas"]):  # type: ignore
            print(
                "%s: M = %s g/mol, m = %s g"
                % (
                    formula,
                    "%.{0}f".format(self.print_precision)
                    % round(self.output["molar masses"][i], self.print_precision),  # type: ignore
                    "%.{0}f".format(self.print_precision)
                    % round(self.output["masses"][i], self.print_precision),  # type: ignore
                )
            )

    def _print_stream(self) -> None:
        """
        Final print stream that can go to different outputs.
        """
        for name, rounded_value in self.rounded_values.items():
            if name == "reaction matrix":
                print(name + ":\n", rounded_value)
            else:
                print(name + ":", rounded_value)
        if self.obj == "ChemicalReaction":
            self._print_additional_reaction_results()

    def print_results(self) -> None:
        """
        Print a final result of calculations in stdout.
        """
        sys.stdout = self.original_stdout
        self._print_stream()

    def write_to_txt(self, filename: str) -> None:
        """
        Export the final result of the calculations in a txt file.

        Arguments:
            filename (str): filename string (should end with .txt)
        """
        if filename == "default":
            filename = self._generate_filename("txt")

        with open(filename, "w", encoding="utf-8") as file:
            sys.stdout = file
            self._print_stream()

        sys.stdout = self.original_stdout

    def dump_to_json(self) -> str:
        """
        Serialization of output into JSON object.

        Returns:
            A JSON-type object
        """
        return json.dumps(self.rounded_values, ensure_ascii=False)

    def write_to_json_file(self, filename: str) -> None:
        """
        Export a final result of calculations in a JSON file.

        Arguments:
            filename (str): filename string (should end with .json)
        """
        if filename == "default":
            filename = self._generate_filename("json")

        with open(filename, "w", encoding="utf-8") as file:
            json.dump(json.loads(self.dump_to_json()), file, ensure_ascii=False)
