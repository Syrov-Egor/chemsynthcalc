import pytest
import csv
import ast
from chemsynthcalc.chemical_reaction import ChemicalReaction
from chemsynthcalc.balancer import Balancer


def test_balancer_wrong_precision():
    reaction = "H2+O2=H2O"
    reaction_obj = ChemicalReaction(reaction)
    with pytest.raises(ValueError):
        Balancer(
            reaction_obj.matrix, len(reaction_obj.decomposed_reaction.reactants), -2
        )


with open("tests/testing_reactions.csv") as csvfile:
    reader = list(csv.reader(csvfile))[1:]
    reactions_set = [(i[0], ast.literal_eval(i[1])) for i in reader]
    inv_set = reactions_set[:100]
    gpinv_set = reactions_set[101:107]
    ppinv_set = reactions_set[108:110]
    comb_set = reactions_set[111:113]


@pytest.mark.parametrize("reaction,coefs", inv_set)
def test_inv_algorithm(reaction: str, coefs: list[int | float]):
    assert ChemicalReaction(reaction).balancer.inv() == coefs


@pytest.mark.parametrize("reaction,coefs", gpinv_set)
def test_gpinv_algorithm(reaction: str, coefs: list[int | float]):
    assert ChemicalReaction(reaction).balancer.gpinv() == coefs


@pytest.mark.parametrize("reaction,coefs", ppinv_set)
def test_ppinv_algorithm(reaction: str, coefs: list[int | float]):
    assert ChemicalReaction(reaction).balancer.ppinv() == coefs


@pytest.mark.parametrize("reaction,coefs", comb_set)
def test_comb_algorithm(reaction: str, coefs: list[int | float]):
    assert ChemicalReaction(reaction).balancer.comb() == coefs
