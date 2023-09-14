import time

# start = time.time()
# time.sleep(1)

from operators import *
from common_eigenstates import *


def unique(list1):
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    return (list(list_set))

############# FUCNTION TO WRITE OUTPUT CSV FILE
def print_to_file(filename, my_tmp, dim):
  f = open(filename, "w")
  f.write("j, i, OP")
  comma=", "
  for i in range(0,dim):
    for j in range(0,dim):
  
      if (abs(my_tmp[i,j]) > 1.e-13) or (i==j):
        print("OP[",i,",", j,"]=", my_tmp[i,j])
        my_row = "\n" + str(j) + comma + str(i) + comma + str(my_tmp[i,j])
        f.write(my_row)
  
  f.close()
#############################

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
#
  Hx = 0.
  Hz = 0.
  dim = 1 << size

# Instanciate Operators
  H_op = hamiltonian(size, J_exch, J_exch, J_exch, Hx, Hz, openBC)
  S_z_op = S_total_z(size)
  T_op = translation(size)

  print("\nIs H real?")
  print(np.all(np.isreal(H_op.matrix)))


  ######## Diogonalize H first ############

  # start = time.time()
  # # Diagonalize Hamiltonian
  # eigenvalue_tmp, eigenvector_tmp = LA.eigh(H_op.matrix, UPLO='U')
  #
  # # reorder and group eigenvalues and eigenvectors
  # H_op_eigenvalue, H_op_eigenvector, H_op_eigenvalue_dict = eigen_dictionary(eigenvalue_tmp, eigenvector_tmp)
  #
  # print("\nrounded H eigenvalues\n")
  # dim_check = 0
  # for ienergy in range(0, len(H_op_eigenvalue)):
  #   E_subset = list(k for k, v in H_op_eigenvalue_dict.items() if v == ienergy)
  #   E_deg = len(E_subset)
  #   print(' %.5f \t %d \t %d' % (H_op_eigenvalue[ienergy], E_deg, ienergy))
  #   dim_check += E_deg
  #
  # print()
  # print(' %s \t %d' % ("Hilbert space dim = ", dim_check))
  #
  # print("\n-------------------------------------\n")
  #
  #
  # H_common_eigenvalue, S_z_common_eigenvalue, common_eigenvector = find_common_eigenstates(S_z_op.matrix, H_op_eigenvalue, H_op_eigenvector, H_op_eigenvalue_dict)
  #
  # end = time.time()
  #
  # print("\nComputational time diagonalizing first H")
  # print(end - start)


  ######## Diogonalize S_z first ############

  start = time.time()

  eigenvalue_tmp =  np.diagonal(S_z_op.matrix)
  eigenvector_tmp = np.diag(np.ones(dim))

  # reorder and group eigenvalues and eigenvectors
  S_z_eigenvalue, S_z_eigenvector, S_z_eigenvalue_dict = eigen_dictionary(eigenvalue_tmp, eigenvector_tmp)

  print("\nrounded S_z eigenvalues\n")
  dim_check = 0
  for ienergy in range(0, len(S_z_eigenvalue)):
    E_subset = list(k for k, v in S_z_eigenvalue_dict.items() if v == ienergy)
    E_deg = len(E_subset)
    print(' %.5f \t %d \t %d' % (S_z_eigenvalue[ienergy], E_deg, ienergy))
    dim_check += E_deg

  print()
  print(' %s \t %d' % ("Hilbert space dim = ", dim_check))

  print("\n-------------------------------------\n")

  S_z_common_eigenvalue, H_common_eigenvalue, common_eigenvector = find_common_eigenstates(H_op.matrix, S_z_eigenvalue, S_z_eigenvector, S_z_eigenvalue_dict)

  end = time.time()

  print("\nComputational time diagonalizing first H")
  print(end - start)


# TODO wirte a test to check that the pairs of eigenvalues obtained diaganalizing H or S_z first are the same
# Reorder by increasing energy
  idx = H_common_eigenvalue.argsort()[::1]

  H_common_eigenvalue_ordered = H_common_eigenvalue[idx]
  S_z_common_eigenvalue_ordered = S_z_common_eigenvalue[idx]
  common_eigenvector_ordered = common_eigenvector[:, idx]


  # ############### ORTH CHECK #############
  # print("orthogonality check")
  # for ialpha in range(0, len(subset)):
  #     for ibeta in range(0, len(subset)):
  #         alpha = subset[ialpha]
  #         beta = subset[ibeta]
  #
  #         test_orth = np.dot(subset_common_eigvector[:, ialpha], subset_common_eigvector[:, ibeta])
  #         if (abs(test_orth) > 1.e-13):
  #             print(alpha, beta, "Op_2 test_orth", test_orth)
  # print()
  # ############### END ORTH CHECK #############



  exit()

# Refactor from here
  output_folder = 'matrix_elements/'
  suffix = 'N=xx.csv'

  E_eigvect = H_op_eigenvector

  OP = T_op  # H_op # S_z_op
  #sp basis
  tmp = np.round(OP.matrix,5)
  #|E> basis
  tmp = np.round(np.matmul(E_eigvect.transpose(), np.matmul(OP.matrix, E_eigvect)),5)

  #|EM> basis
  OP = H_op
  tmp = np.round(np.matmul(S_z_op.E_eigenvector.transpose(), np.matmul(OP.matrix, S_z_op.E_eigenvector)), 5)
  print_to_file(output_folder + "H_EM_basis_" + suffix, tmp, dim)
  OP = S_z_op
  tmp = np.round(np.matmul(S_z_op.E_eigenvector.transpose(), np.matmul(OP.matrix, S_z_op.E_eigenvector)), 5)
  print_to_file(output_folder + "S_z_EM_basis_" + suffix, tmp, dim)
  OP = T_op
  tmp = np.round(np.matmul(S_z_op.E_eigenvector.transpose(), np.matmul(OP.matrix, S_z_op.E_eigenvector)), 5)
  print_to_file(output_folder + "T_EM_basis_" + suffix, tmp, dim)

  exit()
  ########### OUTPUT TERMINAL


  print("S_z on |EM> basis ")
  tmp = np.round(np.matmul(S_z_op.E_eigenvector.transpose(), np.matmul(S_z_op.matrix, S_z_op.E_eigenvector)), 5)
  # print(tmp)
  print("S_z: non-zero matrix elements")
  for i in range(0, dim):
    for j in range(0, dim):
      if (abs(tmp[i, j]) > 1.e-13) or (i == j):
        print("S_z[", i, ",", j, "]=", tmp[i, j])

  print()

  print("H on |EM> basis ")
  tmp = np.round(np.matmul(S_z_op.E_eigenvector.transpose(), np.matmul(H_op.matrix, S_z_op.E_eigenvector)), 5)
  print("H: non-zero matrix elements")
  for i in range(0, dim):
    for j in range(0, dim):
      if (abs(tmp[i, j]) > 1.e-13) or (i == j):
        print("H[", i, ",", j, "]=", tmp[i, j])

  print()
  print("T on |EM> basis ")
  tmp = np.round(np.matmul(S_z_op.E_eigenvector.transpose(), np.matmul(T_op.matrix, S_z_op.E_eigenvector)), 5)
  print("T: non-zero matrix elements")
  for i in range(0, dim):
    for j in range(0, dim):
      if (abs(tmp[i, j]) > 1.e-13) or (i == j):
        print("T[", i, ",", j, "]=", tmp[i, j])


if __name__ == '__main__':
    sys.exit(main())