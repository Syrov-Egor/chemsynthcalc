class MassCalculation():
    '''
    Class for calculation of masses of compounds. 
    Takes reaction coefficients, molar masses list calculated by MolarMassCalculation class,
    index of target compound, target compound mass and rounding_order.
    '''
    def __init__(self, coefficients:list, molar_masses:list, target:int, target_mass:float, rounding_order:int=5):
        self.coefficients:list = coefficients
        self.molar_masses:list = molar_masses
        self.target:int = target
        self.target_mass:float = target_mass
        self.round:int = rounding_order
    
    def calculate_masses(self) -> list:
        '''
        Calculates masses by calculating amount of substance nu (nu=mass/molar mass).
        Coefficients of reaction are normalized to the target. After nu of target compound is
        calculated. it broadcasted to other compound (with respect to its coefficients).
        '''
        if self.coefficients[self.target] != 1.0:
            self.coefficients = [coef/self.coefficients[self.target] for coef in self.coefficients]
        nu = self.target_mass/self.molar_masses[self.target]
        masses = [round(molar*nu*self.coefficients[i], self.round) for i, molar in enumerate(self.molar_masses)]
        return masses
        