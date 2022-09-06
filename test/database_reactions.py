import json
import re
import csv
import sys
import timeit
sys.path.append('./src')
from chemsynthcalc import ChemicalReaction
'''
with open("test\solid-state_dataset_20200713.json", encoding='utf-8') as f:
    data = json.load(f)

data = data.get('reactions')
reactions = [reaction.get('reaction_string') for reaction in data]

with open("test\mined_reactions.csv", 'w', encoding='utf-8') as w:
    writer = csv.writer(w)
    for reaction in reactions:
        writer.writerow([reaction])

in_fnam = "test\mod_mined_reactions.csv"
out_fnam = "test\\allowed_mined_reactions.csv"
with open(in_fnam, newline='', encoding='utf-8') as in_file:
    with open(out_fnam, 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file)
        regex = re.compile(r'[\-x/]')
        search = regex.search
        for row in csv.reader(in_file):
            if not bool(search(row[0])):
                writer.writerow(row)




data = ['4.5 La2O3 + 0.5 Li2CO3 + 6 SiO2 == 1 LiRE9(SiO4)6O2 + 0.5 CO2; RE = La']
data = ['0.5 Li2CO3 + 1.5 MnCO3 + 0.5 O2 + 0.5 TiO2 == 1 LiM0.5Mn1.5O4 + 2 CO2; M = Ti']




def element_substitution(data:list) -> str:
    reaction_list = data[0].split(";")
    if len(reaction_list) > 1:
        if '=' in reaction_list[1]:
            substitute = reaction_list[1]
            atom = re.search(' = (.*)', substitute).group(1)
            sign = re.search(' (.*) = ', substitute).group(1)
            place = re.sub(sign+'(?![a-z])', atom, reaction_list[0])
            place = place.replace("\"", "")
            return place
        else:
            output = reaction_list[0].replace("\"", "")
            return output
    else:
        output = reaction_list[0].replace("\"", "")
        return output 
            

in_fnam = "test\\allowed_mined_reactions.csv"
out_fnam = "test\\final_mined_reactions.csv"
with open(in_fnam, newline='', encoding='utf-8') as in_file:
    with open(out_fnam, 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file)
        data = [row for row in csv.reader(in_file)]
        for row in data:
            new_row = element_substitution(row)
            writer.writerow([new_row])

'''

in_fnam = "test\\final_final_mined_reactions.csv"
with open(in_fnam, encoding='utf-8') as in_file:
    data = [row[0] for row in csv.reader(in_file)]

def calc():
    with open('output.txt', 'w') as file:
        sys.stdout = file
        for reaction in data:
            chemical_reaction = ChemicalReaction(reaction)
            print(chemical_reaction.masses)
