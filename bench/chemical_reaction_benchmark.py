import timeit

from chemsynthcalc import ChemicalReaction


def setup(in_fname: str) -> list[str]:
    with open(in_fname, encoding="utf-8") as reactions:
        data: list[str] = [line.rstrip() for line in reactions]
    return data


def bench(input_list: list[str]):
    with open("bench/reaction_out.txt", "w", encoding="utf-8") as f:
        for reaction in input_list:
            obj = ChemicalReaction(reaction)
            result = str(obj.output_results)
            f.write(reaction + "\n")
            f.write(result + "\n")


input_list = setup("bench/text_mined_reactions.txt")

CYCLES = 1
time_per_cycle = timeit.timeit(lambda: bench(input_list), number=CYCLES) / CYCLES
time_per_reaction = time_per_cycle / len(input_list)

print(f"number of reactions: {len(input_list)}")
print(f"time per cycle: {time_per_cycle} s")
print(f"time per formula: {time_per_reaction * 1000} ms")
