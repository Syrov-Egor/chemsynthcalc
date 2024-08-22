import gc

import numpy as np
import numpy.typing as npt
import scipy #type: ignore


class BalancingAlgorithms:
    def __init__(self, matrix: npt.NDArray[np.float64], separator_pos: int) -> None:
        self.reaction_matrix: npt.NDArray[np.float64] = matrix
        self.reactant_matrix: npt.NDArray[np.float64] = self.reaction_matrix[:, : separator_pos]
        self.product_matrix: npt.NDArray[np.float64] = self.reaction_matrix[:, separator_pos :]

    def _inv_algorithm(self) -> npt.NDArray[np.float64]:
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
        vector = inversed_matrix[:, -zeros_added-1].T
        vector = np.absolute(np.squeeze(np.asarray(vector)))
        vector = vector[vector != 0]
        coefs = np.divide(vector, vector.min())
        return coefs
    
    def _gpinv_algorithm(self) -> npt.NDArray[np.float64]:
        matrix = np.hstack((self.reactant_matrix, -self.product_matrix))
        inverse = scipy.linalg.pinv(matrix)
        a = np.ones((matrix.shape[1], 1))
        i = np.identity(matrix.shape[1])
        coefs = (i - inverse @ matrix) @ a
        return coefs.flat[:]
    
    def _ppinv_algorithm(self) -> npt.NDArray[np.float64]:
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
    
    #!TODO rewrite
    def comb_algorithm(self, max_number_of_iterations: float = 1e8) -> npt.NDArray[np.int32] | None:
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