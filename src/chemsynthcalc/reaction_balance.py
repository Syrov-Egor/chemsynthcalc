import numpy as np
import gc

# import cupy as cp # in case of CuPy combinatorial function
from .chemutils import find_gcd, find_lcm
from fractions import Fraction


class Balancer:
    """
    A class for  balancing chemical equations automatically by different matrix methods.
    Currently implemented: Thorne algorithm (see `inv_algorithm` method for details),
    Risteski general pseudo-inverse algorithm (see `gpinv_algorithm` method for details),
    Risteski partial pseudo-inverse algorithm (see `ppinv_algorithm` method for details),
    and naive combinational search algorithm (see `comb_algorithm` for details).
    Class takes two matrices: matrix of reactants and matrix of products of chemical reaction, derived
    from general reaction matrix. Matrices are in the form of NumPy 2D array.

    Parameters:
    * `reactant_matrix:np.array` - matrix of reactants property generated by `ChemicalReaction` class
    * `product_matrix:np.array` - matrix of products property generated by `ChemicalReaction` class
    * `rounding_order:int` - order of coefficients rounding
    * `intify:bool` - determines whether the coefficients should be integers.
    * `try_comb:bool` - flag, which determines whether an attempt will be made to equalize the reaction
    using the combinatorial method.

    For example:
    ```
    >>>reactant_matrix = ChemicalReaction("H2+O2=H2O").reactant_matrix
    >>>product_matrix = ChemicalReaction("H2+O2=H2O").product_matrix
    >>>Balancer(reactant_matrix, product_matrix, 8, True, False).calculate_coefficients_auto()
    ([2, 1, 2], 'inverse')
    ```
    """

    def __init__(
        self,
        reactant_matrix: np.array,
        product_matrix: np.array,
        rounding_order: int,
        intify: bool,
        try_comb: bool,
        max_comb: int = 1e8,
    ) -> None:
        self.reactant_matrix: np.array = reactant_matrix
        self.product_matrix: np.array = product_matrix
        self.reaction_matrix: np.array = np.hstack(
            (self.reactant_matrix, self.product_matrix)
        )
        self.rounding_order: int = rounding_order
        self.try_comb: bool = try_comb
        self.max_comb: int = max_comb
        self.intify: bool = intify
        self.coef_limit: int = 1000000

    def calculate_coefficients_inv(self) -> list:
        """
        High-level function call to calculate coefficients
        using Thorne algorithm.
        """
        try:
            coefficients = self.inv_algorithm()
            if self.intify:

                intified_coefficients = self.intify_coefficients(
                    coefficients, self.coef_limit
                )
                if self.is_reaction_balanced(
                    self.reactant_matrix, self.product_matrix, intified_coefficients
                ):
                    return intified_coefficients
                elif self.is_reaction_balanced(
                    self.reactant_matrix, self.product_matrix, coefficients
                ):
                    return [int(i) if i.is_integer() else i for i in coefficients]

            else:
                if self.is_reaction_balanced(
                    self.reactant_matrix, self.product_matrix, coefficients
                ):
                    return [int(i) if i.is_integer() else i for i in coefficients]
                else:
                    return None
        except Exception:
            return None

    def calculate_coefficients_gpinv(self) -> list:
        """
        High-level function call to calculate coefficients
        using Risteski algorithm.
        """
        try:
            coefficients = self.gpinv_algorithm()

            if self.intify:
                intified_coefficients = self.intify_coefficients(
                    coefficients, self.coef_limit
                )
                if self.is_reaction_balanced(
                    self.reactant_matrix, self.product_matrix, intified_coefficients
                ):
                    return intified_coefficients
                elif self.is_reaction_balanced(
                    self.reactant_matrix, self.product_matrix, coefficients
                ):
                    return coefficients

            else:
                if self.is_reaction_balanced(
                    self.reactant_matrix, self.product_matrix, coefficients
                ):
                    return coefficients
                else:
                    return None
        except Exception:
            return None

    def calculate_coefficients_ppinv(self) -> list:
        """
        High-level function call to calculate coefficients
        using Risteski algorithm.
        """
        try:
            coefficients = self.ppinv_algorithm()

            if self.intify:
                intified_coefficients = self.intify_coefficients(
                    coefficients, self.coef_limit
                )

                if self.is_reaction_balanced(
                    self.reactant_matrix, self.product_matrix, intified_coefficients
                ):
                    return intified_coefficients
                elif self.is_reaction_balanced(
                    self.reactant_matrix, self.product_matrix, coefficients
                ):
                    return coefficients

            else:
                if self.is_reaction_balanced(
                    self.reactant_matrix, self.product_matrix, coefficients
                ):
                    return coefficients
                else:
                    return None
        except Exception:
            return None

    def calculate_coefficients_comb(self) -> list:
        """
        High-level function call to calculate coefficients
        using combinatorial algorithm.
        """
        try:
            coefficients = self.comb_algorithm(max_number_of_iterations=self.max_comb)
            return coefficients
        except Exception:
            return None

    def calculate_coefficients_auto(self) -> list:
        """
        A method that tries to automatically balance the chemical reaction
        by sequentially applying the three calculation methods to the reaction matrix:
        a Thorne algorithm `inv_algorithm`,
        Risteski algorithm `ppinv_algorithm`, and, if try_comb flag == True,
        a Combinatorial algorithm.
        """
        # try inv
        coefficients = self.calculate_coefficients_inv()
        if coefficients and not any(x <= 0 for x in coefficients):
            return coefficients, "inverse"
        # if failed, try gpinv
        coefficients = self.calculate_coefficients_gpinv()
        if coefficients and not any(x <= 0 for x in coefficients):
            return coefficients, "general pseudoinverse"
        # if failed, try ppinv
        coefficients = self.calculate_coefficients_ppinv()
        if coefficients and not any(x <= 0 for x in coefficients):
            return coefficients, "partial pseudoinverse"
        # if all failed, and try_comb, try comb
        if self.try_comb:
            coefficients = self.calculate_coefficients_comb()
            if coefficients and not any(x <= 0 for x in coefficients):
                return coefficients, "combinatorial"

        print("Cannot automatically balance this reaction")
        return

    @staticmethod
    def is_reaction_balanced(
        reactant_matrix: np.array,
        product_matrix: np.array,
        coefficients: list,
        tolerance: float = 1e-8,
    ) -> bool:
        """
        Checks if reaction is balanced by multiplying reactant matrix and product matrix
        by respective coefficient vector. Method is static to call it outside of balancer
        instance.
        """
        try:
            reactants = np.multiply(
                reactant_matrix.T,
                np.array(coefficients)[: reactant_matrix.shape[1], None],
            )
            products = np.multiply(
                product_matrix.T,
                np.array(coefficients)[reactant_matrix.shape[1] :, None],
            )
            if np.allclose(reactants.sum(axis=0), products.sum(axis=0), rtol=tolerance):
                return True
            else:
                return False
        except Exception:
            return False

    def round_up_coefficients(self, coefficients: list, order: int) -> list:
        """
        Round up the float coefficient if it contains a series
        of the same digit (i.e. 3.00000000009). This function is
        needed to fix inevitable floating point math errors.

        Note:
        There is probably a better and faster way to do it.
        """
        for idx, coefficient in enumerate(coefficients):
            coef_decimal = str(coefficient).split(".")[1]
            coefficients[idx] = coefficient
            for i in range(10):
                if str(i) * (order - 2) in coef_decimal:
                    coefficients[idx] = round(coefficient, order)
                    break
        return coefficients

    def intify_coefficients(self, coefficients: list, limit: int) -> list:
        """
        A function to reduce the coefficients to integers by finding
        greatset common divider.
        """
        initial_coefficients = coefficients
        frac = [Fraction(x).limit_denominator() for x in coefficients]
        vals = [
            int(
                fr.numerator
                * find_lcm([fr.denominator for fr in frac])
                / fr.denominator
            )
            for fr in frac
        ]
        coefficients = [int(val / find_gcd(vals)) for val in vals]
        if any(x > limit for x in coefficients):
            return initial_coefficients
        return coefficients

    def inv_algorithm(self) -> list:
        """
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
        8) Round up float operations errors

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
        """
        reaction_matrix = self.reaction_matrix

        if reaction_matrix.shape[0] > reaction_matrix.shape[1]:
            zero_columns = np.zeros(
                (
                    reaction_matrix.shape[0],
                    reaction_matrix.shape[0] - reaction_matrix.shape[1],
                )
            )
            reaction_matrix = np.hstack((reaction_matrix, zero_columns))

        if reaction_matrix.shape[0] == reaction_matrix.shape[1]:
            p, l, reaction_matrix = np.linalg.svd(reaction_matrix)

        number_of_cols = reaction_matrix.shape[1]
        rank = np.linalg.matrix_rank(reaction_matrix, tol=1e-100)
        nullity = number_of_cols - rank
        augument = np.flip(np.identity(reaction_matrix.shape[1])[:nullity], axis=1)
        augumented_matrix = np.vstack((reaction_matrix, augument))
        if np.where(~augumented_matrix.any(axis=1))[0].size > 0:
            augumented_matrix = augumented_matrix[
                ~np.all(augumented_matrix == 0, axis=1)
            ]
        inversed_matrix = np.linalg.inv(augumented_matrix)
        vector = inversed_matrix[:, -1].T
        vector = np.absolute(np.squeeze(np.asarray(vector)))
        vector = vector[vector != 0]
        coefficients = np.divide(vector, vector.min())
        return self.round_up_coefficients(coefficients.tolist(), self.rounding_order)

    def gpinv_algorithm(self) -> list:
        """
        A reaction matrix pseudoinverse algorithm
        proposed by [Risteski](http://koreascience.or.kr/article/JAKO201314358624990.page).
        There are others articles and methods of chemical
        equation balancing by this author, however, this particular
        algorithm seems to be most convenient for matrix calculations.
        The algorithm can described in steps:

        1) Stack reactant matrix and negative product matrix
        2) Calculate MP pseudoinverse of this matrix
        3) Calculate coefficients by formula:
        x = (I – A+A)a, where x is the coefficients vector,
        I - identity matrix, A+ - MP inverse, A - matrix,
        a - arbitrary vector (in this case, vector of ones).

        Note:
        This method is more general than Thorne method, although it has some
        its own peculiarities. First of all, output of this method is float,
        so, to generate an int coefs list, it needs to be converted, which is
        not always leads to a good result. Secondly, MP pseudoinverse
        is sensetive to row order in the reaction matrix. The rows should
        be ordered by atoms apperances in the reaction string.
        """
        matrix = np.hstack((self.reactant_matrix, -self.product_matrix))
        inverse = np.linalg.pinv(matrix)
        a = np.ones((matrix.shape[1], 1))
        i = np.identity(matrix.shape[1])
        coefs = (i - inverse @ matrix) @ a
        coefs = coefs.flat[:]
        return coefs.tolist()

    def ppinv_algorithm(self) -> list:
        """
        A reaction matrix pseudoinverse algorithm also
        proposed by [Risteski](https://www.koreascience.or.kr/article/JAKO200802727293429.page).
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
        """
        MP_inverse = np.linalg.pinv(self.reactant_matrix)
        g_matrix = (
            np.identity(self.reaction_matrix.shape[0])
            - self.reactant_matrix @ MP_inverse
        )
        g_matrix = g_matrix @ self.product_matrix
        y_multiply = np.linalg.pinv(g_matrix) @ g_matrix
        y_vector = (np.identity(y_multiply.shape[1]) - y_multiply).dot(
            np.ones(y_multiply.shape[1])
        )
        x_multiply = MP_inverse @ self.reactant_matrix
        x_multiply = (
            np.identity(x_multiply.shape[1]) - x_multiply
        ) + MP_inverse @ self.product_matrix @ y_vector.T
        x_vector = x_multiply[0].T
        coefs = np.squeeze(np.asarray(np.hstack((x_vector, y_vector)))).tolist()
        return coefs

    def comb_algorithm(self, max_number_of_iterations: int = 1e8) -> list:
        """
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
        """
        ubyte = 127
        number_of_compounds = self.reaction_matrix.shape[1]
        if number_of_compounds > 10:
            raise ValueError("Sorry, this method is only for n of compound <=10")

        number_of_iterations = int(
            max_number_of_iterations ** (1 / number_of_compounds)
        )

        if number_of_iterations > ubyte:
            number_of_iterations = ubyte

        trans_reaction_matrix = (self.reaction_matrix).T
        lenght = self.reactant_matrix.shape[1]
        old_reactants = trans_reaction_matrix[:lenght].astype("ushort")
        old_products = trans_reaction_matrix[lenght:].astype("ushort")
        for i in range(2, number_of_iterations + 2):
            cart_array = (np.arange(1, i, dtype="ubyte"),) * number_of_compounds
            permuted = np.array(np.meshgrid(*cart_array), dtype="ubyte").T.reshape(
                -1, number_of_compounds
            )
            filter = np.asarray([i - 1], dtype="ubyte")
            permuted = permuted[np.any(permuted == filter, axis=1)]
            # print("calculating max coef %s of %s" % (i-1, number_of_iterations), end='\r', flush=False)
            reactants_vectors = permuted[:, :lenght]
            products_vectors = permuted[:, lenght:]
            del permuted
            reactants = (old_reactants[None, :, :] * reactants_vectors[:, :, None]).sum(
                axis=1
            )
            products = (old_products[None, :, :] * products_vectors[:, :, None]).sum(
                axis=1
            )
            diff = np.subtract(reactants, products)
            del reactants
            del products
            where = np.where(~diff.any(axis=1))[0]
            if np.any(where):
                if where.shape[0] == 1:
                    idx = where
                else:
                    idx = where[0]
                # print("")
                return np.array(
                    np.concatenate(
                        (
                            reactants_vectors[idx].flatten(),
                            products_vectors[idx].flatten(),
                        )
                    )
                ).tolist()
            gc.collect()
        # print("")
        return None

    """
    def comb_algorithm(self, max_number_of_iterations:int = 1e8) -> list:
        '''
        A CuPy GPU-accelerated version of the same combinatorial
        algorithm. Requires [CUDA toolkit] (https://developer.nvidia.com/cuda-toolkit)
        to work. GPU gives around 10x acceleration.

        Finds a solution solution of a Diophantine matrix equation
        by simply enumerating of all possible solutions of number_of_iterations
        coefficients. The solution space is created by Cartesian product
        (in this case, cp.meshgrid function), therefore it is very 
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
        old_reactants = cp.asarray(trans_reaction_matrix[:lenght].astype('ushort'))
        old_products = cp.asarray(trans_reaction_matrix[lenght:].astype('ushort'))
        for i in range(2, number_of_iterations+2):
            cart_array = (cp.arange(1, i, dtype='ubyte'), )*number_of_compounds
            permuted = cp.array(cp.meshgrid(*cart_array), dtype='ubyte').T.reshape(-1,number_of_compounds)
            filter = cp.asarray([i-1], dtype='ubyte')
            permuted = permuted[np.any(permuted==filter, axis=1)]
            print("calculating max coef %s of %s" % (i-1, number_of_iterations), end='\r', flush=True)
            reactants_vectors = cp.asarray(permuted[:, :lenght])
            products_vectors = cp.asarray(permuted[:, lenght:])
            del permuted
            reactants = (old_reactants[None, :, :] * reactants_vectors[:,:, None]).sum(axis=1)
            products = (old_products[None, :, :] * products_vectors[:,:, None]).sum(axis=1)
            diff = cp.subtract(reactants, products)
            del reactants
            del products
            where = cp.where(~diff.any(axis=1))[0]
            if cp.any(where):
                if where.shape[0] == 1:
                    idx = where
                else:
                    idx = where[0]
                print("")
                return cp.array(cp.concatenate((reactants_vectors[idx].flatten(), products_vectors[idx].flatten()))).tolist()
            gc.collect()
        print("")
        print("No solution found")
        return None
    """
