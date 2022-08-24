import numpy as np
import gc
from fractions import Fraction
from .chemutils import find_lcm, find_gcd

class Balancer():
    '''
    A class for  balancing chemical equations automatically by different matrix methods.
    Currently implemented: Thorne algorithm (see `calculate_coefficients_Thorne` method for details),
    Risteski pseudo-inverse algorithm (see `calculate_coefficients_Risteski` method for details),
    and naive combinational search algorithm (see `calculate_coefficients_combinatorial` for details).
    Class takes two matrices: matrix of reactants and matrix of products of chemical reaction, derived
    from general reaction matrix. Matrices are in the form of NumPy 2D array.

    Parameters:
    * `reactant_matrix:np.array`: matrix of reactants property generated by `ChemicalReaction` class
    * `product_matrix:np.array`: matrix of products property generated by `ChemicalReaction` class
    * `rounding_order:int`: order of coefficients rounding
    '''
    def __init__(self, 
    reactant_matrix:np.array, 
    product_matrix:np.array, 
    rounding_order:int) -> None:
        self.reactant_matrix:np.array = reactant_matrix
        self.product_matrix:np.array = product_matrix
        self.reaction_matrix:np.array = np.hstack((self.reactant_matrix, self.product_matrix))
        self.rounding_order:int = rounding_order
        self.coefficient_limit:int = 100000

    def auto_balance_reaction(self) -> list:
        '''
        A method that tries to automatically balance the chemical reaction
        by sequentially applying the three calculation methods to the reaction matrix:
        a Thorne algorithm `calculate_coefficients_Thorne` and
        Risteski algorithm `calculate_coefficients_Risteski`.
        '''
        # try Thorne
        try:
            initial_coefficients = self.calculate_coefficients_Thorne()
            initial_coefficients = self.round_up_coefficients(initial_coefficients, self.rounding_order)
            coefficients = self.intify_coefficients(initial_coefficients, self.coefficient_limit)
            if self.is_reaction_balanced(self.reactant_matrix, self.product_matrix, coefficients):
                return coefficients, "Thorne algorithm"
            elif self.is_reaction_balanced(self.reactant_matrix, self.product_matrix, initial_coefficients):
                return initial_coefficients, "Thorne algorithm"
        except Exception:
            pass
        # if failed, try Risteski
        try:
            initial_coefficients = self.calculate_coefficients_Risteski()
            coefficients = self.intify_coefficients(initial_coefficients, self.coefficient_limit)
            if self.is_reaction_balanced(self.reactant_matrix, self.product_matrix, coefficients):
                return coefficients, "Risteski algorithm"
            elif self.is_reaction_balanced(self.reactant_matrix, self.product_matrix, initial_coefficients):
                return initial_coefficients, "Risteski algorithm"
        except Exception:
            pass
        # if both failed, exit
        print("Can't balance this reaction")
        return

    @staticmethod
    def is_reaction_balanced(reactant_matrix:np.array, product_matrix:np.array, coefficients:list, tolerance:float = 1e-03) -> bool:
        '''
        Checks if reaction is balanced by multiplying reactant matrix and product matrix 
        by respective coefficient vector and. Method is static to call it outside of balancer
        instance.
        '''
        try:
            reactants = np.multiply(reactant_matrix.T, np.array(coefficients)[:reactant_matrix.shape[1], None])
            products = np.multiply(product_matrix.T, np.array(coefficients)[reactant_matrix.shape[1]:, None])
            if np.allclose(reactants.sum(axis=0), products.sum(axis=0), rtol=tolerance):
                return True
            else:
                return False
        except Exception:
            return False

    def round_up_coefficients(self, coefficients:list, order:int) -> list:
        '''
        Round up the float coefficient if it contains a series
        of the same digit (i.e. 3.00000000009). This function is
        needed to fix inevitable floating point math errors.

        Note:
        There is probably a better and faster way to do it.
        '''
        for idx, coefficient in enumerate(coefficients):
            coef_decimal = str(coefficient).split(".")[1]
            coefficients[idx] = coefficient
            for i in range(10):
                if str(i)*(order-2) in coef_decimal:
                    coefficients[idx] = round(coefficient, order)
                    break
        return coefficients

    def intify_coefficients(self, coefficients:list, limit:int) -> list:
        '''
        A function to reduce the coefficients to integers by finding 
        greatset common divider.
        '''
        initial_coefficients = coefficients
        frac = [Fraction(x).limit_denominator() for x in coefficients]
        vals = [int(fr.numerator*find_lcm([fr.denominator for fr in frac])/fr.denominator) for fr in frac]
        coefficients = [int(val/find_gcd(vals)) for val in vals]
        if any(x > limit for x in coefficients):
            return initial_coefficients
        return coefficients
    
    def calculate_coefficients_Thorne(self) -> list:
        '''
        A reaction matrix inverse algorithm proposed by [Thorne](https://arxiv.org/abs/1110.4321).
        The calculation is based on nullity, or dimensionality, of the matrix.
        
        The algorithm can described in steps:
        1) First, check the reaction matrix shape.
        If number of rows in greater then number of columns, add zero
        columns until the matrix becomes square (Note: this is modification
        of original Thorne method described in article).
        2) If reaction matrix is square (which means that number
        of atoms involved is equal to number of compound) than
        we turn matrix in its row-echelon form by singular value
        decomposition.
        3) Calculation of the nullity of the matrix, which is
        basically number of compounds minus rank of the matrix.
        4) Create matrix augumented by nullity number of rows
        of flipped identity matrix. If any rows are zeroes - 
        replace it with identity matrix rows.
        5) Inverse the augumented matrix
        6) Exctract and transpose rightmost column
        7) Normalize this value with absolute min value of vector

        Absolute values of this vector are coefficients for the
        reaction.
        
        Note:
        While this method works great for reactions with 0 and 1
        nullity, it generally cannot work with nullities 2 and higher.
        Thorne claims that for higher nullities, a nullity number
        of vectors should be extracted, and each of them containes
        a set of correct coefficients. However, if number of rows in
        augumentation flipped identity matrix is 2 or more, one can
        easily see that each vector will contain nullity-1 zeroes,
        therefore they cannot be a correct vector of coefficients.
        '''
        reaction_matrix = self.reaction_matrix

        if reaction_matrix.shape[0] > reaction_matrix.shape[1]:
            zero_columns = np.zeros((reaction_matrix.shape[0], reaction_matrix.shape[0]-reaction_matrix.shape[1]))
            reaction_matrix = np.hstack((reaction_matrix,zero_columns))
        
        if reaction_matrix.shape[0] == reaction_matrix.shape[1]:
            p, l, reaction_matrix = np.linalg.svd(reaction_matrix)

        number_of_cols = reaction_matrix.shape[1]
        rank = np.linalg.matrix_rank(reaction_matrix, tol = 1e-100)
        nullity = number_of_cols - rank
        augument = np.flip(np.identity(reaction_matrix.shape[1])[:nullity], axis = 1)
        augumented_matrix = np.vstack((reaction_matrix, augument))
        if np.where(~augumented_matrix.any(axis=1))[0]:
            augumented_matrix = augumented_matrix[~np.all(augumented_matrix == 0, axis=1)]
        inversed_matrix = np.linalg.inv(augumented_matrix)
        vector = inversed_matrix[:, -1].T
        vector = np.absolute(np.squeeze(np.asarray(vector)))
        vector = vector[vector != 0]
        coefficients = np.divide(vector, vector.min())
        return coefficients.tolist()

    def calculate_coefficients_Risteski(self) -> list:
        '''
        A reaction matrix pseudoinverse algorithm proposed by [Risteski](https://www.koreascience.or.kr/article/JAKO200802727293429.page).
        There are others articles and methods of chemical
        equation balancing by this author, however, this particular 
        algorithm seems to be most convenient for matrix calculations.
        The method is founded on virtue of the solution of a 
        Diophantine matrix equation by using of a Moore-Penrose 
        pseudoinverse matrix.

        The algorithm can described in steps:
        1) Take the Moore-Penrose pseudoinverse of reactant matrix
        2) Create a G matrix in the form of (I-AA^-)B, where
        I is the identity matrix, A is the reactant matrix, A^- is
        the MP pseudoinverse of A and B is the product matrix.
        3) Then, vector y (coefficients of products) is equal to
        (I-G^-G)u. 
        4) Vector x (coefficients of reactants) is equal to 
        A^-By + (I-A^-A)v, where u and v are columns of ones.

        Note:
        This method is more general than Thorne method, although it has some
        its own peculiarities. First of all, output of this method is float,
        so, to generate an int coefs list, it needs to be converted, which is
        not always leads to a good result. Secondly, MP pseudoinverse
        is sensetive to row order in the reaction matrix. The rows should
        be ordered by atoms apperances in the reaction string.
        '''
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

    def calculate_coefficients_combinatorial(self, max_number_of_iterations:int = 1e8) -> list:
        '''
        Finds a solution solution of a Diophantine matrix equation
        by simply enumerating of all possible solutions of number_of_iterations
        coefficients. The solution space is created by Cartesian product
        (in this case, np.meshgrid function), therefore it is very 
        limited by memory. There must a better, clever and fast solution 
        to this.

        Note:
        All possible variations of coefficients vectors are
        `combinations = max_coefficients**number_of_compounds`
        therefore this method is most effective for reaction with
        small numbers of compounds.
        '''
        ubyte = 127
        number_of_compounds = self.reaction_matrix.shape[1]
        if number_of_compounds>10:
            raise ValueError("Sorry, this method is for n of compound <=10")

        number_of_iterations = int(max_number_of_iterations**(1/number_of_compounds))

        if number_of_iterations > ubyte:
            number_of_iterations = ubyte
    
        trans_reaction_matrix = (self.reaction_matrix).T
        lenght = self.reactant_matrix.shape[1]
        old_reactants = trans_reaction_matrix[:lenght].astype('ushort')
        old_products = trans_reaction_matrix[lenght:].astype('ushort')
        for i in range(2, number_of_iterations+2):
            cart_array = (np.arange(1, i, dtype='ubyte'), )*number_of_compounds
            permuted = np.array(np.meshgrid(*cart_array), dtype='ubyte').T.reshape(-1,number_of_compounds)
            filter = np.asarray([i-1], dtype='ubyte')
            permuted = permuted[np.in1d(permuted[:, 1], filter)]
            print("calculating %s of %s row" % (i-1, number_of_iterations), end='\r', flush=True)
            reactants_vectors = permuted[:, :lenght]
            products_vectors = permuted[:, lenght:]
            del permuted
            reactants = (old_reactants[None, :, :] * reactants_vectors[:,:, None]).sum(axis=1)
            products = (old_products[None, :, :] * products_vectors[:,:, None]).sum(axis=1)
            diff = np.subtract(reactants, products)
            print(diff.dtype)
            del reactants
            del products
            where = np.where(~diff.any(axis=1))[0]
            if np.any(where):
                if where.shape[0] == 1:
                    idx = where
                else:
                    idx = where[0]
                print("")
                return np.array(np.concatenate((reactants_vectors[idx].flatten(), products_vectors[idx].flatten()))).tolist(), "Combinatorial algorithm"
            gc.collect()
        print("")
        print("No solution found")
        return None, "Combinatorial algorithm"

    def calculate_coefficients_combinatorial_old(self, max_number_of_iterations:int = 1e8) -> list:
        '''
        Finds a solution solution of a Diophantine matrix equation
        by simply enumerating of all possible solutions of number_of_iterations
        coefficients. The solution space is created by Cartesian product
        (in this case, np.meshgrid function), therefore it is very 
        limited by memory. There must a better, clever and fast solution 
        to this.

        Note:
        All possible variations of coefficients vectors are
        `combinations = max_coefficients**number_of_compounds`
        therefore this method is most effective for reaction with
        small numbers of compounds.
        '''

        ubyte = 127
        number_of_compounds = self.reaction_matrix.shape[1]
        if number_of_compounds>10:
            raise ValueError("Sorry, this method is for n of compounds <= 10")

        number_of_iterations = int(max_number_of_iterations**(1/number_of_compounds))

        if number_of_iterations > ubyte:
            number_of_iterations = ubyte
    
        trans_reaction_matrix = (self.reaction_matrix).T.astype('short')
        lenght = self.reactant_matrix.shape[1]
        for i in range(2, number_of_iterations+2):
            cart_array = (np.arange(1, i, dtype='ubyte'), )*number_of_compounds
            permuted = np.array(np.meshgrid(*cart_array), dtype='ubyte').T.reshape(-1,number_of_compounds)
            filter = np.asarray([i-1], dtype='ubyte')
            vectors = permuted[np.in1d(permuted[:, 1], filter)]
            print("calculating %s of %s max coef" % (i-1, number_of_iterations), end='\r', flush=True)
            for i in range(vectors.shape[0]):
                summ = trans_reaction_matrix * vectors[i][:, None]
                if np.array_equal(summ[:lenght].sum(axis=0), summ[lenght:].sum(axis=0)):
                    print("")
                    return vectors[i].tolist(), "Combinatorial algorithm"
        print("")
        print("No solution found")
        return None, "Combinatorial algorithm"