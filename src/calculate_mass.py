class MassCalculation():
    def __init__(self, coefficients:list, molar_masses:list, target:int, target_mass:float):
        self.coefficients = coefficients
        self.molar_masses = molar_masses
        self.target = target
        self.target_mass = target_mass
        self.round = 5
    
    def calculate_masses(self) -> list:
        if self.coefficients[self.target] != 1.0:
            self.coefficients = [coef/self.coefficients[self.target] for coef in self.coefficients]
        nu = self.target_mass/self.molar_masses[self.target]
        masses = [round(molar*nu*self.coefficients[i], self.round) for i, molar in enumerate(self.molar_masses)]
        return masses
        