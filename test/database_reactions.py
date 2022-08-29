import json
import csv
import sys
sys.path.append('./src')
from chemsynthcalc import ChemicalReaction
'''
with open("tests\solid-state_dataset_2019-06-27_upd.json", encoding='utf-8') as f:
    data = json.load(f)

reactions = [point.get('reaction_string') for point in data]

with open("tests\mined_reactions.csv", 'w', encoding='utf-8') as w:
    writer = csv.writer(w)
    for reaction in reactions:
        writer.writerow([reaction])
'''
out_fnam = "tests\good_mined_reactions.csv"
with open("tests\mod_mined_reactions.csv", newline='') as in_file:
    reader = csv.reader(in_file)
    reactions = [row for row in reader]

allowed = []
for reaction in reactions:

    try:
        react_obj = ChemicalReaction(reaction[0])
        allowed.append(react_obj)
    except Exception:
        pass

print(len(allowed))