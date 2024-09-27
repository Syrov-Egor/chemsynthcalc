import sys
import time
import json

import numpy as np

from .utils import round_dict_content


class ChemicalOutput:
    def __init__(
        self, output: dict[str, object], print_precision: int, obj: str
    ) -> None:
        if print_precision > 0:
            self.print_precision: int = print_precision
        else:
            raise ValueError("precision <= 0")

        if obj == "formula" or obj == "reaction":
            self.obj = obj
        else:
            raise ValueError(f"No such object: {obj}")

        self.output: dict[str, object] = output
        self.rounded_values: dict[str, object] = self._round_values()
        self.original_stdout = sys.stdout

    def _round_values(self) -> dict[str, object]:
        rounded_dict: dict[str, object] = {}
        for name, value in self.output.items():
            if isinstance(value, float):
                rounded_value = round(value, self.print_precision)
            elif isinstance(value, dict):
                rounded_value = round_dict_content(value, self.print_precision)  # type: ignore
            elif name == "masses":
                rounded_value = [round(v, self.print_precision) for v in value]  # type: ignore
            elif name == "reaction matrix":
                rounded_value = np.array2string(value) # type: ignore
            else:
                rounded_value = value

            rounded_dict.update({name: rounded_value})  # type: ignore

        return rounded_dict

    def _generate_filename(self, file_type: str) -> str:
        if self.obj == "formula":
            filename: str = (
                f"CSC_{self.obj}_{self.output.get("formula")}_{time.time_ns()}.{file_type}"
            )
        else:
            filename: str = (
                f"CSC_{self.obj}_{self.output.get("target")}_{time.time_ns()}.{file_type}"
            )

        return filename

    def _print_additional_reaction_results(self) -> None:
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
        for name, rounded_value in self.rounded_values.items():
            if name == "reaction matrix":
                print(name + ":\n", rounded_value)
            else:
                print(name + ":", rounded_value)
        if self.obj == "reaction":
            self._print_additional_reaction_results()

    def print_results(self) -> None:
        sys.stdout = self.original_stdout
        self._print_stream()

    def write_to_txt(self, filename: str) -> None:
        if filename == "default":
            filename = self._generate_filename("txt")

        with open(filename, "w") as file:
            sys.stdout = file
            self._print_stream()

        sys.stdout = self.original_stdout

    def dump_to_json(self) -> str:
        return json.dumps(self.rounded_values, ensure_ascii=False)

    def write_to_json_file(self, filename: str) -> None:
        if filename == "default":
            filename = self._generate_filename("json")

        with open(filename, "w", encoding="utf-8") as file:
            json.dump(json.loads(self.dump_to_json()), file, ensure_ascii=False)
