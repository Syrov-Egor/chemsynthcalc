import numpy as np
from scipy.linalg import lu
from fractions import Fraction
from .chemutils import find_lcm, find_gcd

class Balancer():
    '''
    A class for automatically balancing chemical equations by different matrix methods.
    Currently implemented Thorne algorithm (see calculate_coefficients_Thorne method for details),
    Risteski pseudo-inverse algorithm (see calculate_coefficients_Risteski method for details),
    and naive combinational search algorithm (see calculate_coefficients_combinatorial for details).
    Class takes two matrices: matrix of reactants and matrix of products of chemical reaction, derived
    from general reaction matrix. Matrices are in the form of NumPy 2D array.
    '''
    def __init__(self, reactant_matrix:np.array, product_matrix:np.array, rounding_order:int = 5) -> None:
        self.reactant_matrix:np.array = reactant_matrix
        self.product_matrix:np.array = product_matrix
        self.reaction_matrix:np.array = np.hstack((self.reactant_matrix, self.product_matrix))
        self.rounding_order:int = rounding_order

    def auto_balance_reaction(self) -> list:
        '''
        A method that tries to automatically balance the chemical reaction
        by sequentially applying the three calculation methods to the reaction matrix:
        a Thorne algorithm, a Risteski algorithm and combinatorial search algorithm.
        '''
        try:
            coefficients = self.calculate_coefficients_Thorne()
            coefficients = self.roundUpCoefficients(coefficients, self.rounding_order)
            if self.is_reaction_balanced(self.reactant_matrix, self.product_matrix, coefficients):
                return [int(i) if i.is_integer() else i for i in coefficients], "Thorne algorithm"
        except:
            try:
                coefficients = self.calculate_coefficients_Risteski()
                #coefficients = self.roundUpCoefficients(coefficients, self.rounding_order)
                frac = [Fraction(x).limit_denominator() for x in coefficients]
                vals = [int(fr.numerator*find_lcm([fr.denominator for fr in frac])/fr.denominator) for fr in frac]
                coefficients = [int(val/find_gcd(vals)) for val in vals]
                if self.is_reaction_balanced(self.reactant_matrix, self.product_matrix, coefficients):
                    return coefficients, "Risteski algorithm"
            except:
                try:
                    coefficients = self.calculate_coefficients_combinatorial()
                    if self.is_reaction_balanced(self.reactant_matrix, self.product_matrix, coefficients):
                        return coefficients, "Combinational search"
                except:
                    print("Cannot equalize this reaction...")
                    return []
    
    @staticmethod
    def is_reaction_balanced(reactant_matrix:np.array, product_matrix:np.array, coefficients:list, tolerance:float = 1e-03) -> bool:
        '''

        '''
        reactants = np.multiply(reactant_matrix.T, np.array(coefficients)[:reactant_matrix.shape[1], None])
        products = np.multiply(product_matrix.T, np.array(coefficients)[reactant_matrix.shape[1]:, None])
        if np.allclose(reactants.sum(axis=0), products.sum(axis=0), rtol=tolerance):
            return True
        else:
            return False

    def roundUpCoefficients(self, coefficients:list, order:int) -> list:
        for idx, coefficient in enumerate(coefficients):
            coef_decimal = str(coefficient).split(".")[1]
            coefficients[idx] = coefficient
            for i in range(10):
                if str(i)*(order) in coef_decimal:
                    coefficients[idx] = round(coefficient, order)
                    break
        return coefficients
    
    def calculate_coefficients_Thorne(self) -> list:
        reaction_matrix = self.reaction_matrix
        #print("initial matrix", reaction_matrix)
        if reaction_matrix.shape[0] > reaction_matrix.shape[1]:
            zero_columns = np.zeros((reaction_matrix.shape[0], reaction_matrix.shape[0]-reaction_matrix.shape[1]))
            reaction_matrix = np.hstack((reaction_matrix,zero_columns))
        
        if reaction_matrix.shape[0] == reaction_matrix.shape[1]:
            p, l, reaction_matrix = lu(reaction_matrix)

        number_of_cols = reaction_matrix.shape[1]
        #print("reaction matrix", reaction_matrix)
        rank = np.linalg.matrix_rank(reaction_matrix, tol = 1e-100)
        #print("rank of matrix", rank)
        nullity = number_of_cols - rank
        #print("nullity", nullity)
        augument = np.flip(np.identity(reaction_matrix.shape[1])[:nullity], axis = 1)
        augumented_matrix = np.vstack((reaction_matrix, augument))
        if np.where(~augumented_matrix.any(axis=1))[0]:
            augumented_matrix = augumented_matrix[~np.all(augumented_matrix == 0, axis=1)]
        #print("augumented_matrix", augumented_matrix)
        inversed_matrix = np.linalg.inv(augumented_matrix)
        #print("inverted matrix", inversed_matrix)
        vector = inversed_matrix[:, -1].T
        vector = np.absolute(np.squeeze(np.asarray(vector)))
        vector = vector[vector != 0]
        #print("vector", vector)
        coefficients = np.divide(vector, vector.min())
        return coefficients.tolist()

    def calculate_coefficients_Risteski(self):
        MP_inverse = np.linalg.pinv(self.reactant_matrix)
        g_matrix = np.identity(self.reaction_matrix.shape[0])-self.reactant_matrix@MP_inverse
        g_matrix = g_matrix@self.product_matrix
        y_multiply = np.linalg.pinv(g_matrix)@g_matrix
        y_vector = (np.identity(y_multiply.shape[1])-y_multiply).dot(np.ones(y_multiply.shape[1]))
        x_multiply = MP_inverse@self.reactant_matrix
        x_multiply = (np.identity(x_multiply.shape[1])-x_multiply) + MP_inverse@self.product_matrix@y_vector.T
        x_vector = x_multiply[0].T
        coefs = np.squeeze(np.asarray(np.hstack((x_vector, y_vector)))).tolist()
        return coefs

    def calculate_coefficients_combinatorial(self, number_of_iterations:int = 10) -> list:
        trans_reaction_matrix = (self.reaction_matrix).T
        lenght = self.reactant_matrix.shape[1]
        number_of_compounds = self.reaction_matrix.shape[1]
        old_reactants = trans_reaction_matrix[:lenght]
        old_products = trans_reaction_matrix[lenght:]
        for i in range(2, number_of_iterations+2):
            cart_array = (np.arange(1, i),)*number_of_compounds
            permuted = np.array(np.meshgrid(*cart_array)).T.reshape(-1,number_of_compounds)
            print("calculating %s of %s row..." % (i-1, number_of_iterations))
            for item in permuted:
                if i-1 in item:
                    with np.errstate(all='ignore'):
                        reactants = np.multiply(old_reactants, np.array(item)[:lenght, None])
                        products = np.multiply(old_products, np.array(item)[lenght:, None])
                        if np.array_equal(reactants.sum(axis=0), products.sum(axis=0)):
                            return np.array(item).tolist()
        print("No solution")
        return None