import gc

import numpy as np
import numpy.typing as npt
import scipy  # type: ignore


class BalancingAlgorithms:
    """
    A collection of functions for balancing chemical reactions

    Currently implemented: Thorne algorithm (see
    [_inv_algorithm][chemsynthcalc.balancing_algos.BalancingAlgorithms._inv_algorithm] method for details),
    Risteski general pseudo-inverse algorithm (see
    [_gpinv_algorithm][chemsynthcalc.balancing_algos.BalancingAlgorithms._gpinv_algorithm] method for details),
    Risteski partial pseudo-inverse algorithm (see
    [_ppinv_algorithm][chemsynthcalc.balancing_algos.BalancingAlgorithms._ppinv_algorithm] method for details),
    and naive combinational search algorithm (see
    [_comb_algorithm][chemsynthcalc.balancing_algos.BalancingAlgorithms._comb_algorithm] method for details).

    Parameters:
        matrix (npt.NDArray[np.float64]): Reaction matrix
        separator_pos (int): Position of the reaction separator (usually the separator is "=")

    Attributes:
        reactant_matrix (npt.NDArray[np.float64]): A matrix of the left part of the equation
        product_matrix (npt.NDArray[np.float64]): A matrix of the right part of the equation

    Note:
        Why use
        [scipy.linalg.pinv](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.pinv.html),
        when
        [numpy.linalg.pinv](https://numpy.org/doc/stable/reference/generated/numpy.linalg.pinv.html)
        is doing the same thing and does not require the whole SciPy import?

        There are some peculiar reaction cases where
        (especially for [_ppinv_algorithm][chemsynthcalc.balancing_algos.BalancingAlgorithms._ppinv_algorithm] method)
        the results for [numpy.linalg.pinv](https://numpy.org/doc/stable/reference/generated/numpy.linalg.pinv.html)
        differs from system to system (np version, OS, python version etc.). My understanding is that the cause of
        this behaviour lies in small differences for pinv algorithm in numpy C-libraries and BLAS-libraries,
        hence the difference.
        To avoid this, a more consistent method
        [scipy.linalg.pinv](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.pinv.html) was used.
    """

    def __init__(self, matrix: npt.NDArray[np.float64], separator_pos: int) -> None:
        self.separator_pos = separator_pos
        self.reaction_matrix: npt.NDArray[np.float64] = matrix
        self.reactant_matrix: npt.NDArray[np.float64] = self.reaction_matrix[
            :, : self.separator_pos
        ]
        self.product_matrix: npt.NDArray[np.float64] = self.reaction_matrix[
            :, self.separator_pos :
        ]

    def _inv_algorithm(self) -> npt.NDArray[np.float64]:
        """Matrix inverse algorithm for reaction balancing.

        A reaction matrix inverse algorithm proposed by [Thorne](https://arxiv.org/abs/1110.4321).
        The calculation is based on the nullity, or dimensionality, of the matrix.

        The algorithm can be described in steps:

        1) If the number of rows is greater than the number of columns, \
        add zero columns until the matrix becomes square \
        (Note: this is a modification of the original \
        Thorne method described in the article).

        2) If reaction matrix is square (which means that the number \
        of atoms involved is equal to the number of compounds) than \
        we turn matrix in its row-echelon form by singular value \
        decomposition.

        3) Calculation of the nullity of the matrix, which is \
        basically number of compounds minus rank of the matrix.

        4) Create a  matrix augumented by nullity number of rows \
        of flipped identity matrix. If any rows are zeros, \
        replace them with identity matrix rows.

        5) Inverse the augumented matrix.

        6) Exctract and transpose rightmost column.

        7) Normalize this value with the absolute min value of the vector.

        8) Round up float operations errors.

        The absolute values of this vector are coefficients of the
        reaction.

        Note:
            While this method works great for reactions with 0 and 1
            nullity, it generally cannot work with nullities 2 and higher.
            Thorne claims that for higher nullities, a nullity number
            of vectors should be extracted, and each of them contains
            a set of correct coefficients. However, if number of rows in
            the flipped augmentation identity matrix is 2 or more, one can
            easily see that each vector will contain nullity-1 zeroes,
            therefore they cannot be a correct vector of coefficients.

        Returns:
            A 1D NumPy array of calculated coefficients
        """
        reaction_matrix = self.reaction_matrix

        if reaction_matrix.shape[0] > reaction_matrix.shape[1]:
            zeros_added = reaction_matrix.shape[0] - reaction_matrix.shape[1]
            zero_columns = np.zeros(
                (
                    reaction_matrix.shape[0],
                    zeros_added,
                )
            )
            reaction_matrix = np.hstack((reaction_matrix, zero_columns))
        else:
            zeros_added = 0

        if reaction_matrix.shape[0] == reaction_matrix.shape[1]:
            _, _, reaction_matrix = np.linalg.svd(reaction_matrix)

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
        vector = inversed_matrix[:, -zeros_added - 1].T
        vector = np.absolute(np.squeeze(np.asarray(vector)))
        vector = vector[vector != 0]
        coefs = np.divide(vector, vector.min())
        return coefs

    def _gpinv_algorithm(self) -> npt.NDArray[np.float64]:
        """Matrix gerenal pseudoinverse algorithm for reaction balancing.

        A reaction matrix pseudoinverse algorithm
        proposed by [Risteski](http://koreascience.or.kr/article/JAKO201314358624990.page).
        There are other articles and methods of chemical
        equation balancing by this author, however, this particular
        algorithm seems to be most convenient for matrix calculations.
        The algorithm can be described in steps:

        1) Stack reactant matrix and negative product matrix.

        2) Calculate MP pseudoinverse of this matrix.

        3) Calculate coefficients by formula:
        x = (I – A+A)a, where x is the coefficients vector,
        I - identity matrix, A+ - MP inverse, A - matrix,
        a - arbitrary vector (in this case, vector of ones).

        Note:
            This method is more general than Thorne's method, although it has some
            peculiarities of its own. First of all, the output of this method is float array,
            so, to generate an int coefs list, it needs to be converted, which is
            not always leads to a good result. Secondly, MP pseudoinverse
            is sensetive to row order in the reaction matrix. The rows should
            be ordered by atoms apperances in the reaction string.

        Returns:
            A 1D NumPy array of calculated coefficients
        """
        matrix = np.hstack((self.reactant_matrix, -self.product_matrix))
        inverse = scipy.linalg.pinv(matrix)
        a = np.ones((matrix.shape[1], 1))
        i = np.identity(matrix.shape[1])
        coefs = (i - inverse @ matrix) @ a
        return coefs.flat[:]

    def _ppinv_algorithm(self) -> npt.NDArray[np.float64]:
        """
        Matrix partial pseudoinverse algorithm for reaction balancing.

        A reaction matrix pseudoinverse algorithm also
        proposed by [Risteski](https://www.koreascience.or.kr/article/JAKO200802727293429.page).
        The method is founded on virtue of the solution of a
        Diophantine matrix equation by using of a Moore-Penrose
        pseudoinverse matrix.

        The algorithm can be described in steps:

        1) Take the Moore-Penrose pseudoinverse of the reactant matrix.

        2) Create a G matrix in the form of (I-AA^-)B, where
        I is the identity matrix, A is the reactant matrix, A^- is
        the MP pseudoinverse of A and B is the product matrix.

        3) Then, the vector y (coefficients of products) is equal to
        (I-G^-G)u.

        4) Vector x (coefficients of reactants) is equal to
        A^-By + (I-A^-A)v, where u and v are columns of ones.

        Note:
            While this algorithm and
            [_gpinv_algorithm][chemsynthcalc.balancing_algos.BalancingAlgorithms._gpinv_algorithm]
            are very similar, there are some differences in output results.
            This method exists mostly for legacy purposes, like balancing
            some reactions according to [Risteski](https://www.koreascience.or.kr/article/JAKO200802727293429.page).

        Returns:
            A 1D NumPy array of calculated coefficients
        """
        MP_inverse = scipy.linalg.pinv(self.reactant_matrix)
        g_matrix = (
            np.identity(self.reaction_matrix.shape[0])
            - self.reactant_matrix @ MP_inverse
        )
        g_matrix = g_matrix @ self.product_matrix
        y_multiply = scipy.linalg.pinv(g_matrix) @ g_matrix
        y_vector = (np.identity(y_multiply.shape[1]) - y_multiply).dot(
            np.ones(y_multiply.shape[1])
        )
        x_multiply = MP_inverse @ self.reactant_matrix
        x_multiply = (
            np.identity(x_multiply.shape[1]) - x_multiply
        ) + MP_inverse @ self.product_matrix @ y_vector.T
        x_vector = x_multiply[0].T
        coefs = np.squeeze(np.asarray(np.hstack((x_vector, y_vector))))
        return coefs

    def _comb_algorithm(
        self, max_number_of_iterations: float = 1e8
    ) -> npt.NDArray[np.int32] | None:
        """
        Matrix combinatorial algorithm for reaction balancing.

        Finds a solution solution of a Diophantine matrix equation
        by simply enumerating of all possible solutions of number_of_iterations
        coefficients. The solution space is created by Cartesian product
        (in this case, *np.meshgrid* function), and therefore it is very
        limited by memory. There must a better, clever and fast solution
        to this!

        Important:
            Only for integer coefficients less than 128. Only for reactions
            with total compound count <=10.
            A GPU-accelerated version of this method can be done by importing
            CuPy and replacing np. with cp.

        Note:
            All possible variations of coefficients vectors are
            combinations = max_coefficients**number_of_compounds,
            therefore this method is most effective for reaction with
            small numbers of compounds.

        Returns:
            A 1D NumPy array of calculated coefficients of None if can't compute
        """
        byte = 127
        number_of_compounds = self.reaction_matrix.shape[1]
        if number_of_compounds > 10:
            raise ValueError("Sorry, this method is only for n of compound <=10")

        number_of_iterations = int(
            max_number_of_iterations ** (1 / number_of_compounds)
        )

        if number_of_iterations > byte:
            number_of_iterations = byte

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
                )
            gc.collect()
        # print("")
        return None
