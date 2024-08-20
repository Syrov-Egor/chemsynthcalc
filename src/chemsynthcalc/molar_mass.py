from .periodic_table import PeriodicTable
from .formula_parser import ChemicalFormulaParser


class MolarMassCalculation:
    def __init__(self, parsed_formula: dict[str, float]) -> None:
        self.parsed_formula: dict[str, float] = parsed_formula
        self.p_table = PeriodicTable().p_table

    def _calculate_atomic_masses(self) -> list[float]:
        masses: list[float] = []
        for atom, weight in self.parsed_formula.items():
            p_table_row: tuple[float, float, str] | None = self.p_table.get(atom)
            if p_table_row is not None:
                atom_mass: float = p_table_row[0]
            else:
                raise ValueError
            masses.append(atom_mass * weight)
        return masses

    def calculate_molar_mass(self) -> float:
        return sum(self._calculate_atomic_masses())

    def calculate_mass_percent(self) -> dict[str, float]:
        atomic_masses: list[float] = self._calculate_atomic_masses()
        molar_mass: float = self.calculate_molar_mass()
        percents: list[float] = [atomic / molar_mass * 100 for atomic in atomic_masses]
        return dict(zip(self.parsed_formula.keys(), percents))

    def calculate_atomic_percent(self) -> dict[str, float]:
        values: list[float] = list(self.parsed_formula.values())
        atomic: list[float] = [value / sum(values) * 100 for value in values]
        return dict(zip(self.parsed_formula.keys(), atomic))

    def calculate_oxide_percent(self) -> dict[str, float]:
        mass_percents: list[float] = list(self.calculate_mass_percent().values())
        oxides: list[tuple[str, str, float]] = []
        for i, atom in enumerate(self.parsed_formula.keys()):
            if atom != "O":
                p_table_row: tuple[float, float, str] | None = self.p_table.get(atom)
                if p_table_row is not None:
                    oxide_label: str = p_table_row[2]
                    oxides.append((atom, oxide_label, mass_percents[i]))

        oxide_percents: list[float] = []

        for oxide in oxides:
            parsed_oxide: dict[str, float] = ChemicalFormulaParser(
                oxide[1]
            ).parse_formula()
            p_table_row: tuple[float, float, str] | None = self.p_table.get(oxide[0])
            oxide_mass: float = MolarMassCalculation(
                parsed_oxide
            ).calculate_molar_mass()
            atomic_oxide_coef: float | None = parsed_oxide.get(oxide[0])
            if p_table_row is not None and atomic_oxide_coef is not None:
                atomic_mass: float = p_table_row[0]
                conversion_factor: float = oxide_mass / atomic_mass / atomic_oxide_coef
                oxide_percents.append(oxide[2] * conversion_factor)

        normalized_oxide_percents: list[float] = [
            x / sum(oxide_percents) * 100 for x in oxide_percents
        ]
        oxide_labels: list[str] = [oxide[1] for oxide in oxides]

        return dict(zip(oxide_labels, normalized_oxide_percents))
