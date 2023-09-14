# USAGE
# $ pytest test_operators.py -v
#
from operators import *
from common_eigenstates import *
import numpy as np

size = 6
J_exch = 0.25
Hx = 0.
Hz = 0.
openBC = False # open boundary conditions set to True should break est_commutation_H_T()
H_op = hamiltonian(size, J_exch, J_exch, J_exch, Hx, Hz, openBC)
S_x_op = S_total_x(size)
S_y_op = S_total_y(size)
S_z_op = S_total_z(size)
T_op = translation(size)

def test_H_herimtean():
    np.testing.assert_array_equal(H_op.matrix, H_op.matrix.transpose())

def test_commutation_H_S_x():
#   Operators H and S_z commute
    P1 = np.dot(S_x_op.matrix, H_op.matrix)
    P2 = np.dot(H_op.matrix, S_x_op.matrix)
    np.testing.assert_array_equal(P1, P2)

def test_commutation_H_S_y():
#   Operators H and S_y commute
    P1 = np.dot(S_y_op.matrix, H_op.matrix)
    P2 = np.dot(H_op.matrix, S_y_op.matrix)
    np.testing.assert_array_equal(P1, P2)

def test_commutation_H_S_z():
#   Operators H and S_z commute
    P1 = np.dot(H_op.matrix, S_z_op.matrix)
    P2 = np.dot(S_z_op.matrix, H_op.matrix)
    np.testing.assert_array_equal(P1, P2)

def test_commutation_H_T():
#   Operators H and T commute
    # TODO check why it fails for size = 8 with Hz = 0.2
    P1 = np.dot(H_op.matrix, T_op.matrix)
    P2 = np.dot(T_op.matrix, H_op.matrix)
    np.testing.assert_array_equal(P1, P2)
    #np.testing.assert_allclose(P1, P2, rtol=1e-5, atol=0) # otherwise it may fail because of Max relative difference: 1.44431353e-16

def test_commutation_S_x_T():
#   Operators S_x and T commute
    P1 = np.dot(S_x_op.matrix, T_op.matrix)
    P2 = np.dot(T_op.matrix, S_x_op.matrix)
    np.testing.assert_array_equal(P1, P2)

def test_commutation_S_y_T():
#   Operators S_y and T commute
    P1 = np.dot(S_y_op.matrix, T_op.matrix)
    P2 = np.dot(T_op.matrix, S_y_op.matrix)
    np.testing.assert_array_equal(P1, P2)