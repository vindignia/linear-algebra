import sys
import numpy as np
from numpy import linalg as LA
from mlxtend.math import vectorspace_orthonormalization


def are_same(A, B):
    # Check if A and B are square and have the same dimension. If so, check if they are identical and returns a boolean
    if (np.shape(A)[0] != np.shape(A)[1]):
        sys.exit('A is not a square matrix')
    if (np.shape(B)[0] != np.shape(B)[1]):
        sys.exit('B is not a square matrix')
    if (np.shape(A)[0] != np.shape(B)[0]):
        sys.exit('A and B have not the same size')

    # Convert the lists to NumPy arrays
    A = np.array(A)
    B = np.array(B)

    # Check if the arrays are equal element-wise
    return np.array_equal(A, B)


def do_commute(A, B):
    # Check if the product of two square matrices commute
    C = np.dot(A, B)
    D = np.dot(B, A)
    return are_same(C, D)


def order_eigenvalue_vectors(eigenvalue, eigenvector, asc=True):
    if not asc:
        idx = eigenvalue.argsort()[::-1]
    else:
        idx = eigenvalue.argsort()[::1]

    eigenvalue_ordered = eigenvalue[idx]
    eigenvector_ordered = eigenvector[:, idx]

    return eigenvalue_ordered, eigenvector_ordered


def eigen_dictionary(eigenvalue, eigenvector, asc=True):
    # order eigenvector and eigenvalues
    eigenvalue_tmp, eigenvector_tmp = order_eigenvalue_vectors(eigenvalue, eigenvector, asc)

    eigenvalue_round = np.round(eigenvalue_tmp.real, 5)
    eigenvalue, indices, degeneracy = np.unique(eigenvalue_round, return_counts=True, return_index=True)

    eigenvalue_dict = {}
    ii = 0
    for ieig_1 in range(0, len(eigenvalue)):
        # normalize eigenvectors associated with degenerate eigenvalues
        deg_subset = range(ii, ii + degeneracy[ieig_1])  # it should work also when deg = 1
        subset_eigvect = vectorspace_orthonormalization(eigenvector_tmp[:, deg_subset])

        for i in range(0, degeneracy[ieig_1]):
            eigenvalue_dict[ii] = ieig_1
            eigenvector[:, ii] = subset_eigvect[:, i]
            ii += 1

    return eigenvalue, eigenvector, eigenvalue_dict


def find_common_eigenstates(B, A_eigenvalue, A_eigenvector, A_eigenvalue_dict):
    common_eigenvector = np.zeros(np.shape(B))
    A_common_eigenvalue = np.zeros(shape=len(A_eigenvalue_dict))
    B_common_eigenvalue = np.zeros(shape=len(A_eigenvalue_dict))

    ii = 0  # label of eigevectors associated to pairs of common eigenvalues

    for ieig_1 in range(0, len(A_eigenvalue)):

        subset = list(k for k, v in A_eigenvalue_dict.items() if v == ieig_1)

        subset_A_eigenvector = A_eigenvector[:, subset]
        subset_common_eigvector_not_norm = np.zeros(shape=(len(subset), len(subset)))
        subset_common_eigvector = np.zeros(shape=(len(subset), len(subset)))

        # matrix elements of B on the H_op eigenstates of the present subset
        B_subset_tmp = np.matmul(B, subset_A_eigenvector)
        B_subset = np.matmul(subset_A_eigenvector.transpose(), B_subset_tmp)

        if len(subset) > 1:
            # diagonalize <beta|B|alpha>
            B_subset_eigvalue_not_norm, subset_common_eigvector_not_norm = LA.eigh(B_subset)
            # round and group
            B_subset_eigvalue, indices, subset_B_degeneracy = np.unique(np.round(B_subset_eigvalue_not_norm, 5),
                                                                        return_counts=True,
                                                                        return_index=True)
        else:  # energy singlet
            B_subset_eigvalue, subset_B_degeneracy = np.round(B_subset, 5), 1

        jj = 0
        for ieig_2 in range(0, len(B_subset_eigvalue)):  # distinct B eigenvalues in the current subset

            if (len(subset) > 1) and (subset_B_degeneracy[ieig_2] > 1):
                # normalize eigenvectors associated with degenerate eigenvalues
                B_deg_subset = range(jj, jj + subset_B_degeneracy[ieig_2])
                subset_common_eigvector[:, B_deg_subset] = vectorspace_orthonormalization(
                    subset_common_eigvector_not_norm[:, B_deg_subset])
                # eigenvalues
                eig_2_deg_subset = range(ii + jj, ii + jj + subset_B_degeneracy[ieig_2])
                B_common_eigenvalue[eig_2_deg_subset] = B_subset_eigvalue[ieig_2]
                for i in range(0, subset_B_degeneracy[ieig_2]):
                    jj += 1

            elif (len(subset) > 1 and (subset_B_degeneracy[ieig_2] == 1)):
                # take care of the non-deg M_eigenstates NOT energy singlet!!!
                subset_common_eigvector[:, jj] = subset_common_eigvector_not_norm[:, jj]
                B_common_eigenvalue[ii + jj] = B_subset_eigvalue[ieig_2]
                jj += 1

            else:  # singlet
                # In this case np assign a shape (1,1) to B_subset_eigvalue which has to be accessed to be considered a scalar
                B_common_eigenvalue[ii] = B_subset_eigvalue[0, 0]

        if len(subset) > 1:
            subset_common_tmp = np.matmul(subset_A_eigenvector, subset_common_eigvector)
            subset_common = range(ii, ii + len(subset))
            common_eigenvector[:, subset_common] = subset_common_tmp
            A_common_eigenvalue[subset_common] = A_eigenvalue[ieig_1]
            ii += len(subset)
        else:  # singlet
            common_eigenvector[:, ii] = A_eigenvector[:, ii]
            A_common_eigenvalue[ii] = A_eigenvalue[ieig_1]
            ii += 1

    return A_common_eigenvalue, B_common_eigenvalue, common_eigenvector
