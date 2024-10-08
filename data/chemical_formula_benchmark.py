import csv
import timeit

from chemsynthcalc import ChemicalFormula
from chemsynthcalc import ChemicalReaction
from chemsynthcalc.chemutils import stripe_formula_from_coefficients

def setup(in_fname: str) -> list[str]:
    with open(in_fname, encoding='utf-8') as reactions:
        reader = csv.reader(reactions)
        data: list[str] = [row[0] for row in reader]

    all_fomulas: list[str] = []
    for reaction in data:
        decomposed = ChemicalReaction(reaction).compounds
        decomposed = [stripe_formula_from_coefficients(formula)[1] for formula in decomposed]
        all_fomulas.extend(decomposed)

    return list(set(all_fomulas))

def bench(input_list: list[str]):
    for formula in input_list:
        obj = ChemicalFormula(formula)
        result = obj.output_results
        print(result)


       

input_list = setup("data/text_mined_reactions.csv")

with open("out.txt", "w", encoding="utf-8") as f:
    for line in input_list:
        f.write(line + "\n")

CYCLES = 10
time_per_cycle = timeit.timeit(lambda:bench(input_list), number=CYCLES) / CYCLES
time_per_formula = time_per_cycle / len(input_list)

print(f"number of formulas: {len(input_list)}")
print(f"time per cycle: {time_per_cycle} s")
print(f"time per formula: {time_per_formula * 1000} ms")