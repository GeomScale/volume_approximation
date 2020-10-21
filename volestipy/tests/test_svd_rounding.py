#!/usr/bin/python3.6

import numpy as np
import volestipy
from volestipy import HPolytope

if __name__ == "__main__":
   dim = 5
   A = np.zeros((2*dim, dim), dtype=np.float)
   A[0:dim] = np.eye(dim)
   A[dim:] -=  np.eye(dim,dim, dtype=np.float)
   b = np.ones(2*dim, dtype=np.float)

   hp = HPolytope(A,b)

   print("The polytope:")
   print("Dimensions: {}".format(hp.dimensions))
   print("A matrix:")
   print(hp.A)
   print("b vector:")
   print(hp.b)

   rounding_returns = ["new_A","new_b","T_matrix","shift","round_val"]

   svd_rounding = hp.rounding_svd()

   print("\n\n ** This is the rounding svd output\n")
   for i in range(len(rounding_returns)):
      print("\n" + rounding_returns[i] + ":")
      print(svd_rounding[0])

