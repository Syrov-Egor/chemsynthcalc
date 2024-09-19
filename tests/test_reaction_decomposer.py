from chemsynthcalc.reaction_decomposer import ReactionDecomposer


def test_split_coefficient_from_formula():
    reaction = "2H2+O2=2H2O"
    assert ReactionDecomposer(reaction).initial_coefficients == [2.0, 1.0, 2.0]
