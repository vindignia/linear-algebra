from operators import *
from common_eigenstates import *
import pandas as pd


def compute_matrix_elements(M, basis_vect):
  return np.round(np.matmul(basis_vect.transpose(), np.matmul(M, basis_vect)), 5)


def print_to_file(filename, my_tmp, dim):
  f = open(filename, "w")
  f.write("j, i, OP")
  comma=", "
  for i in range(0,dim):
    for j in range(0,dim):
  
      if abs(my_tmp[i,j]) > 1.e-13 or (i==j):
        my_row = "\n" + str(j) + comma + str(i) + comma + str(my_tmp[i,j])
        f.write(my_row)
  
  f.close()


def main():

  if (len(sys.argv) > 1):
    size = int(sys.argv[1])  # number of spins 1/2 in the ring
  else:
    size = 4

  if (len(sys.argv) > 2):
    openBC = sys.argv[2]
  else:
    openBC = False

  if (len(sys.argv) > 3):
    eps = float(sys.argv[3])
  else:
    eps = 0.

  J_exch = 0.25  # Exchange interaction [Kelvin] time spin modulus J/4
  Hx = 0.
  Hz = 0.

# Instanciate Operators
  H_op = hamiltonian(size, J_exch, J_exch, J_exch, Hx, Hz, openBC)
  S_z_op = S_total_z(size)
  T_op = translation(size)

  # Diagonalize Hamiltonian
  eigenvalue_tmp, eigenvector_tmp = LA.eigh(H_op.matrix, UPLO='U')

  # reorder and group eigenvalues and eigenvectors
  H_op_eigenvalue, H_op_eigenvector, H_op_eigenvalue_dict = eigen_dictionary(eigenvalue_tmp, eigenvector_tmp)

  print("\nH eigenvalues\tmultiplicity\n")
  dim_check = 0
  for ienergy in range(0, len(H_op_eigenvalue)):
    H_subset = list(k for k, v in H_op_eigenvalue_dict.items() if v == ienergy)
    H_deg = len(H_subset)
    print(' %.5f \t %d' % (H_op_eigenvalue[ienergy], H_deg))
    dim_check += H_deg

  print()
  print(' %s \t %d' % ("Hilbert space dim =", dim_check))

  print("\n-------------------------------------")

  H_common_eigenvalue, S_z_common_eigenvalue, common_eigenvector = find_common_eigenstates(S_z_op.matrix, H_op_eigenvalue, H_op_eigenvector, H_op_eigenvalue_dict)

  df = pd.DataFrame({'H_eig':H_common_eigenvalue, 'S_z_eig':S_z_common_eigenvalue})
  df = df.sort_values(['H_eig','S_z_eig'], ascending=[True, True])

  print("\nH and S_z common eigenvalues\n")
  print(df)

  output_folder = 'matrix_elements/'
  dim = 1 << size
  suffix = 'N'+  str(dim) + '.csv'

  print("\nMatrix elements on the common eigenvector basis written in the folder " + output_folder +"\n")

  # matrix elements on the common |H, S_z> basis
  tmp = compute_matrix_elements(H_op.matrix, common_eigenvector)
  print_to_file(output_folder + "H_comm_basis_" + suffix, tmp, dim)

  tmp = compute_matrix_elements(S_z_op.matrix, common_eigenvector)
  print_to_file(output_folder + "S_z_comm_basis_" + suffix, tmp, dim)

  tmp = compute_matrix_elements(T_op.matrix, common_eigenvector)
  print_to_file(output_folder + "T_comm_basis_" + suffix, tmp, dim)

if __name__ == '__main__':
    sys.exit(main())