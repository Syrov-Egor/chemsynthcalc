import pytest

from chemsynthcalc.molar_mass import MolarMassCalculation

molar_mass_test_data: list[tuple[dict[str, float], float]] = [
    ({"H": 2.0, "O": 1.0}, 18.015),
    ({"N": 2.0, "H": 10.0, "S": 1.0, "O": 5.0}, 150.149),
    ({"K": 1.2, "Na": 0.8, "S": 1.0, "O": 4.0}, 161.365415424),
    ({"K": 4.0, "Mg": 2.0, "S": 6.0, "O": 24.0, "Ho": 2.0}, 1111.198658),
]


@pytest.mark.parametrize("parsed_formula,molar_mass", molar_mass_test_data)
def test_molar_mass(parsed_formula: dict[str, float], molar_mass: float):
    assert MolarMassCalculation(parsed_formula).calculate_molar_mass() == molar_mass


mass_percent_test_data: list[tuple[dict[str, float], dict[str, float]]] = [
    ({"H": 2.0, "O": 1.0}, {"H": 11.19067443796836, "O": 88.80932556203163}),
    (
        {"N": 2.0, "H": 10.0, "S": 1.0, "O": 5.0},
        {
            "N": 18.657466916196576,
            "H": 6.713331424118708,
            "S": 21.35212355726645,
            "O": 53.277078102418265,
        },
    ),
    (
        {"K": 1.2, "Na": 0.8, "S": 1.0, "O": 4.0},
        {
            "K": 29.075375213902195,
            "Na": 11.397619109196414,
            "S": 19.86794996670129,
            "O": 39.659055710200114,
        },
    ),
    (
        {"K": 4.0, "Mg": 2.0, "S": 6.0, "O": 24.0, "Ho": 2.0},
        {
            "K": 14.074171065098604,
            "Mg": 4.374555319162381,
            "S": 17.31103602538728,
            "O": 34.5551173262846,
            "Ho": 29.685120264067127,
        },
    ),
]


@pytest.mark.parametrize("parsed_formula,mass_percent", mass_percent_test_data)
def test_mass_percent(parsed_formula: dict[str, float], mass_percent: dict[str, float]):
    assert MolarMassCalculation(parsed_formula).calculate_mass_percent() == mass_percent


atomic_percent_test_data: list[tuple[dict[str, float], dict[str, float]]] = [
    ({"H": 2.0, "O": 1.0}, {"H": 66.66666666666666, "O": 33.33333333333333}),
    (
        {"N": 2.0, "H": 10.0, "S": 1.0, "O": 5.0},
        {
            "N": 11.11111111111111,
            "H": 55.55555555555556,
            "S": 5.555555555555555,
            "O": 27.77777777777778,
        },
    ),
    (
        {"K": 1.2, "Na": 0.8, "S": 1.0, "O": 4.0},
        {
            "K": 17.142857142857142,
            "Na": 11.428571428571429,
            "S": 14.285714285714285,
            "O": 57.14285714285714,
        },
    ),
    (
        {"K": 4.0, "Mg": 2.0, "S": 6.0, "O": 24.0, "Ho": 2.0},
        {
            "K": 10.526315789473683,
            "Mg": 5.263157894736842,
            "S": 15.789473684210526,
            "O": 63.1578947368421,
            "Ho": 5.263157894736842,
        },
    ),
]


@pytest.mark.parametrize("parsed_formula,atomic_percent", atomic_percent_test_data)
def test_atomic_percent(
    parsed_formula: dict[str, float], atomic_percent: dict[str, float]
):
    assert (
        MolarMassCalculation(parsed_formula).calculate_atomic_percent()
        == atomic_percent
    )


oxide_percent_test_data: list[tuple[dict[str, float], dict[str, float]]] = [
    ({"H": 2.0, "O": 1.0}, {"H2O": 100.0}),
    (
        {"N": 2.0, "H": 10.0, "S": 1.0, "O": 5.0},
        {"NO2": 35.09929732740271, "H2O": 34.36114777487011, "SO3": 30.53955489772719},
    ),
    (
        {"K": 1.2, "Na": 0.8, "S": 1.0, "O": 4.0},
        {
            "K2O": 35.02423357043221,
            "Na2O": 15.363524680216429,
            "SO3": 49.61224174935137,
        },
    ),
    (
        {"K": 4.0, "Mg": 2.0, "S": 6.0, "O": 24.0, "Ho": 2.0},
        {
            "K2O": 16.713129118300564,
            "MgO": 7.151185901417123,
            "SO3": 42.61382168343717,
            "Ho2O3": 33.52186329684514,
        },
    ),
]


@pytest.mark.parametrize("parsed_formula,oxide_percent", oxide_percent_test_data)
def test_oxide_percent(
    parsed_formula: dict[str, float], oxide_percent: dict[str, float]
):
    assert (
        MolarMassCalculation(parsed_formula).calculate_oxide_percent() == oxide_percent
    )
