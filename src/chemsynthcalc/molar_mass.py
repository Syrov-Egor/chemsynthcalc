from .periodic_table import PeriodicTable

class MolarMassCalculation:
    def __init__(self, parsed_formula: dict[str, float]) -> None:
        self.parsed_formula: dict[str, float] = parsed_formula

    def _calculate_atomic_masses(self) -> list[float]:
        masses: list[float] = []
        for atom, weight in self.parsed_formula.items():
            p_table_row: tuple[float, float, str] | None = PeriodicTable().p_table.get(atom)
            if p_table_row is not None:
                atom_mass: float = p_table_row[0]
            else:
                raise ValueError
            masses.append(atom_mass*weight)
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