import sys
import numpy as np
import time
from numpy import linalg as LA

start = time.time()
time.sleep(1)

from operators import hamiltonian
from common_eigenstates import  *

def array_spin_eigenvalues(size, J_exch, openBC):

  Hx = 0.
  Hz = 0.
  H_op = hamiltonian(size,J_exch,J_exch,J_exch,Hx,Hz,openBC)

  # Diagonalize Hamiltonian
  eigenvalue_tmp, eigenvector_tmp = LA.eigh(H_op.matrix, UPLO='U')

  # reorder and group eigenvalues and eigenvectors
  H_op_eigenvalue, H_op_eigenvector, H_op_eigenvalue_dict = eigen_dictionary(eigenvalue_tmp, eigenvector_tmp)

  print("\nrounded H eigenvalues\n")
  dim_check = 0
  for ienergy in range(0, len(H_op_eigenvalue)):
    E_subset = list(k for k, v in H_op_eigenvalue_dict.items() if v == ienergy)
    E_deg = len(E_subset)
    print(' %.5f \t %d \t %d' % (H_op_eigenvalue[ienergy], E_deg, ienergy))
    dim_check += E_deg

  print()
  print(' %s \t %d' % ("Hilbert space dim = ", dim_check))


def main():


  if (len(sys.argv) > 1):
    size = int(sys.argv[1])  # number of spins 1/2 in the ring
  else:
    size = 3

  if (len(sys.argv) > 2):
    openBC = sys.argv[2]
  else:
    openBC = False

  if (len(sys.argv) > 3):
    eps = float(sys.argv[3])
  else:
    eps = 0.

  J_exch = 0.25  # Exchange interaction [Kelvin] time spin modulus J/4

  array_spin_eigenvalues(size, J_exch, openBC)


if __name__ == '__main__':
    sys.exit(main())