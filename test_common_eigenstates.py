# USAGE
# $ pytest test_common_eigenstates.py -v
#
from operators import *
from common_eigenstates import *
import numpy as np

size = 6
J_exch = 0.25
Hx = 0.
Hz = 0.
H_op = hamiltonian(size, J_exch, J_exch, J_exch, Hx, Hz, False)
S_z_op = S_total_z(size)
S_x_op = S_total_z(size)
T_op = translation(size)


def find_H_eigenvalues_eigenvectors(H_op_matrix):
    eigenvalue_tmp, eigenvector_tmp = LA.eigh(H_op_matrix, UPLO='U')
    return eigen_dictionary(eigenvalue_tmp, eigenvector_tmp)


def find_S_z_eigenvalues_eigenvectors(size):
    S_z_op = S_total_z(size)
    eigenvalue_tmp = np.diagonal(S_z_op.matrix)
    dim = 1 << size
    eigenvector_tmp = np.diag(np.ones(dim))
    return eigen_dictionary(eigenvalue_tmp, eigenvector_tmp)


def test_consistency_eigenvalues_eigenvectors():
    H_op_eigenvalue, H_op_eigenvector, H_op_eigenvalue_dict = find_H_eigenvalues_eigenvectors(H_op.matrix)

    H_common_eigenvalue, S_x_common_eigenvalue, common_eigenvector = find_common_eigenstates(S_x_op.matrix,
                                                                                             H_op_eigenvalue,
                                                                                             H_op_eigenvector,
                                                                                             H_op_eigenvalue_dict)

    test_S_x_eig = np.array([])
    test_H_eig = np.array([])

    for ii in range(0, len(H_common_eigenvalue)):
        vec = common_eigenvector[:, ii]
        tmp_S_x_eig = np.round(np.matmul(vec.transpose(), np.matmul(S_x_op.matrix, vec)), 5)
        tmp_H_eig = np.round(np.matmul(vec.transpose(), np.matmul(H_op.matrix, vec)), 5)
        test_S_x_eig = np.append(test_S_x_eig, tmp_S_x_eig)
        test_H_eig = np.append(test_H_eig, tmp_H_eig)

    np.testing.assert_array_equal(H_common_eigenvalue, test_H_eig, err_msg='eigenvalues inconsistent with eigenvectors',
                                  verbose=True)
