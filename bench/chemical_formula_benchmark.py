import csv
import timeit

from chemsynthcalc import ChemicalFormula
from chemsynthcalc.reaction_decomposer import ReactionDecomposer


def setup(in_fname: str) -> list[str]:
    with open(in_fname, encoding="utf-8") as reactions:
        reader = csv.reader(reactions)
        data: list[str] = [row[0] for row in reader]

    all_fomulas: list[str] = []
    for reaction in data:
        decomposed = ReactionDecomposer(reaction).compounds
        all_fomulas.extend(decomposed)

    return list(set(all_fomulas))


def bench(input_list: list[str]):
    with open("bench/formula_out.txt", "w", encoding="utf-8") as f:
        for formula in input_list:
            obj = ChemicalFormula(formula)
            result = str(obj.output_results)
            f.write(formula + "\n")
            f.write(result + "\n")


input_list = setup("data/text_mined_reactions.csv")

CYCLES = 1
time_per_cycle = timeit.timeit(lambda: bench(input_list), number=CYCLES) / CYCLES
time_per_formula = time_per_cycle / len(input_list)

print(f"number of formulas: {len(input_list)}")
print(f"time per cycle: {time_per_cycle} s")
print(f"time per formula: {time_per_formula * 1000} ms")
