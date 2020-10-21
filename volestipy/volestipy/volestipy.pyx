#!python
#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False

# Global dependencies
import os
import sys
import numpy as np
cimport numpy as np
from libcpp cimport bool
from cpython cimport bool

# For the preprocess step, we need the following dependencies
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

# For the read the json format BIGG files function
import json
import scipy.io

# For the scaling routine
from scipy.sparse import diags
import math

# ----------------------------------------------------------------------------------

# Set the time
def get_time_seed():
   import random
   import time
   return int(time.time())

## Read a Bigg file and get the necessary A and b
# The .json format case
def read_json_file(input_file):

   with open(input_file, 'r') as f:

      distros_dict = json.load(f)

      reactions_list = distros_dict['reactions']

      metabolites = []
      reactions = []

      for reaction in reactions_list:

         metabolites_dic = reaction['metabolites']
         reaction_name = reaction['id']
         reactions.append(reaction_name)

         for metabolite in metabolites_dic.keys():
            if metabolite not in metabolites:
               metabolites.append(metabolite)

      list_of_reaction_lists = []
      vector_of_ubs = []
      vector_of_lbs = []

      for reaction in reactions_list:

         ub = float(reaction['upper_bound']) ; vector_of_ubs.append(ub)
         lb = float(reaction['lower_bound']) ; vector_of_lbs.append(lb)

         metabolites_dic = reaction['metabolites']
         reaction_vector = []

         for term in range(len(metabolites)):

            metabolite = metabolites[term]

            if metabolite in metabolites_dic.keys():

               reaction_vector.append(metabolites_dic[metabolite])
            else:
               reaction_vector.append(0)

         list_of_reaction_lists.append(reaction_vector)

   f.close()

   # Build function's output; first the A matrix
   n = len(list_of_reaction_lists)
   A = np.zeros((2*n, n), dtype=np.float)
   A[0:n] = np.eye(n)
   A[n:] -=  np.eye(n,n, dtype=np.float)

   # Now, the b vector
   vector_of_lbs = [-x for x in vector_of_lbs]
   b = np.asarray(vector_of_ubs + vector_of_lbs)

   # The Aeq matrix
   Aeq = np.asarray(list_of_reaction_lists)
   Aeq = np.transpose(Aeq)

   # And the beq vector
   m = len(metabolites)
   beq = np.zeros(m)
   
   # Make everything C contigeous
   A = np.asarray(A, dtype = 'float')
   A = np.ascontiguousarray(A, dtype='float')
   b = np.asarray(b, dtype = 'float')
   b = np.ascontiguousarray(b, dtype='float')
   Aeq = np.asarray(Aeq, dtype = 'float')
   Aeq = np.ascontiguousarray(Aeq, dtype='float')
   beq = np.asarray(beq, dtype = 'float')
   beq = np.ascontiguousarray(beq, dtype='float')

   return A, b, Aeq, beq, metabolites, reactions

# The .mat format case
def read_mat_file(input_file):

   data_from_mat = scipy.io.loadmat(input_file)

   species_name = ''
   for key in data_from_mat.keys():
      if key[0] != "_":
         species_name = key

   species = data_from_mat[species_name]
   list_of_lists = species.tolist()

   counter = 0

   metabolites = []

   for element in list_of_lists[0][0]:

      if counter == 0:

         m =len(element)

         for i in element:

            metabolite = i[0][0]

            if metabolite not in metabolites:
               metabolites.append(metabolite)

      if counter == 7:
         reactions_list = element.tolist()
         reactions = [reaction[0][0] for reaction in reactions_list]
      if counter == 11:
         lb_tmp = element
         lb_tmp = lb_tmp.tolist()
      if counter == 12:
         ub_tmp = element
      if counter == 10:
         Aeq = element

      counter += 1

   # Build function's output; first the A matrix
   n = len(ub_tmp)
   A = np.zeros((2*n, n), dtype=np.float)
   A[0:n] = np.eye(n)
   A[n:] -=  np.eye(n,n, dtype=np.float)

   # Now, the b vector
   ub = [i[0] for i in ub_tmp]
   lb = [-x[0] for x in lb_tmp]
   b = np.asarray(ub + lb)

   # The Aeq matrix
   Aeq = np.asarray(Aeq)

   # And the beq vector
   beq = np.zeros(m)

   return A, b, Aeq, beq, metabolites, reactions

# Build a Python functionto pre-process the metabolic network; meaning to remove really "small" facets.
# This function will be implemented by making use of the Gurobi solver
def pre_process(A, b, Aeq, beq):

   d = A.shape[1] ; m = Aeq.shape[0] ; n = Aeq.shape[1]
   Aeq_new = Aeq
   beq_new = beq
   A_new = np.zeros((0,d))
   b_new = []    # this need to be a vector; we do not know its length
   
   min_fluxes = []
   max_fluxes = []

   try:

      # To avoid printint the output of the optimize() function of Gurobi, we need to set an environment like this
      with gp.Env(empty=True) as env:
          env.setParam('OutputFlag', 0)
          env.start()

          with gp.Model(env=env) as model:

            # Create variables
            x = model.addMVar(shape = d, vtype = GRB.CONTINUOUS , name = "x", lb = -GRB.INFINITY, ub = GRB.INFINITY)

            # Make sparse Aeq
            # Aeq = np.array(Aeq, dtype=float)
            # Aeq_sparse = sp.csr_matrix(Aeq.astype(np.float))
            Aeq_sparse = sp.csr_matrix(Aeq)

            # Make A sparse
            # A = np.array(A, dtype=float)
            # A_sparse = sp.csr_matrix(A.astype(np.float))
            A_sparse = sp.csr_matrix(A)            

            # Set the b and beq vectors as numpy vectors
            b = np.array(b)
            beq = np.array(beq)

            # Add constraints
            model.addMConstrs(Aeq_sparse, x, '=', beq, name = "c")

            # Update the model to include the constraints added
            model.update()

            #######################
            # After getting the constraints you need to add the bounds; ObjBound might work:
            # https://www.gurobi.com/documentation/9.0/refman/objbound.html#attr:ObjBound
            # to start with, avoid ObjBound and do that the same way as Aeq but with unequalities this time
            #######################

            # Add constraints for the uneqalities of A
            model.addMConstrs(A_sparse, x, '<', b, name = "d")

            # Update the model with the extra constraints and then print it
            model.update()
            model.display()

            # Loop through the lines of the A matrix, set objective function for each and run the model
            for i in range(A.shape[0]):
            
               if Aeq_new.shape[0] < Aeq_new.shape[1] - 10:

                  # Set the ith row of the A matrix as the objective function
                  objective_function = A[i,]
   
                  # Set the objective function in the model
                  model.setMObjective(None, objective_function, 0.0, None, None, x, GRB.MINIMIZE)
                  model.update()
   
                  # Optimize model
                  model.optimize ()
   
                  # If optimized
                  status = model.status
                  if status == GRB.OPTIMAL:
   
                     # Get the max objective value
                     max_objective = model.getObjective().getValue()
                     max_fluxes.append(max_objective)
   
                  # Likewise, for the minimum
                  objective_function = np.asarray([-x for x in objective_function])
                  model.setMObjective(None, objective_function, 0.0, None, None, x, GRB.MINIMIZE)
                  model.update()
                  model.optimize()
   
                  # Again if optimized
                  status = model.status
                  if status == GRB.OPTIMAL:
   
                     # Get the max objective value
                     min_objective = model.getObjective().getValue()
                     min_fluxes.append(min_objective)
   
                  # Calculate the width
                  # width = abs(max_objective + min_objective) / np.linalg.norm(A[i,])
                  width = abs(max_objective - min_objective)                  
   
                  # Check whether we keep or not the equality
                  if width < 1e-07:
                     Aeq_new = np.vstack((Aeq_new, A[i,]))
                     # beq_new = np.append(beq_new, max_objective)
                     beq_new = np.append(beq_new, min(max_objective, min_objective))
   
                  else:
                     A_new = np.vstack((A_new, A[i,]))
                     b_new = np.append(b_new, b[i])
                     
               else:
                  A_new = np.vstack((A_new, A[i,]))
                  b_new = np.append(b_new, b[i])                  
               

            # The np.vstack() creates issues on changing contiguous c orded of np arrays; here we fix this
            Aeq_new = np.ascontiguousarray(Aeq_new, dtype=np.dtype)
            A_new = np.ascontiguousarray(A_new, dtype=np.dtype)

            # Furthremore, we need to have float64 in all numpy arrays
            Aeq_new = Aeq_new.astype('float64')
            A_new = A_new.astype('float64')

            # Make lists of fluxes numpy arrays
            min_fluxes = np.asarray(min_fluxes) ; max_fluxes = np.asarray(max_fluxes)

            # And now we export the pre-processed elements on .npy files to load and use them anytime
            np.save('A_preprocessed.npy', A_new) ; np.save('b_preprocessed.npy', b_new)
            np.save('Aeq_preprocessed.npy', Aeq_new) ; np.save('beq_preprocessed.npy', beq_new)
            np.save('min_fluxes.npy', min_fluxes), np.save('max_fluxes.npy', max_fluxes)
            
            # Return a tupple including the new A, b, Aeq and beq
            return A_new, b_new, Aeq_new, beq_new, min_fluxes, max_fluxes

   # Print error messages
   except gp . GurobiError as e :
      print ("Error code " + str( e . errno ) + ": " + str( e ))
   except AttributeError :
      print ("Encountered an attribute error ")

# This is a function to get the maximum ball included in the full dimensional polytope
def get_max_ball(A_full_dim, b_full_dim):

   extra_column = []

   m = A_full_dim.shape[0]
   n = A_full_dim.shape[1]

   for i in range(A_full_dim.shape[0]):
      entry = np.linalg.norm(A_full_dim[i,])
      extra_column.append(entry)

   column = np.asarray(extra_column)
   A_expand = np.c_[A_full_dim, column]

   with gp.Env(empty=True) as env:
      env.setParam('OutputFlag', 0)
      env.start()

      d = A_expand.shape[1]

      with gp.Model(env=env) as model:

         # Create variable
         x = model.addMVar(shape = d, vtype = GRB.CONTINUOUS , name = "x", lb = -GRB.INFINITY, ub = GRB.INFINITY)
         model.update()

         # Make A_full_dim sparse
         # A_expand = np.array(A_expand, dtype=float)
         A_expand_sparse = sp.csr_matrix(A_expand.astype(np.float))

         # Add constraints
         model.addMConstrs(A_expand_sparse, x, '<', b_full_dim, name = "c")
         model.update()

         # Set the ith row of the A matrix as the objective function
         a = np.ones((1,n+1))
         b = np.zeros((1,n))
         a[:,:-1] = b
         objective_function = a[0]

         # Set the objective function in the model
         model.setMObjective(None, objective_function, 0.0, None, None, x, GRB.MAXIMIZE)
         model.update()

         # Optimize model
         model.optimize ()

         # Get the solution returned
         vars = model.getVars()

         # Get the center point and the radius of max ball from the solution of LP; its last element is the radius
         point = []
         for i in range(len(vars)):
            if i == len(vars) - 1:
               r = vars[i].x
            else:
               value = vars[i].x
               point.append(value)

         # And check whether its value is negative
         if r < 0 :
            print ("The radius calculated has negative value. The polytope is infeasible or something went wrong with the solver")
         else:
            return point, r

# Map the points samples on the (rounded) full dimensional polytope, back to the initial one
def map_samples_on_initial_polytope(samples, T, T_shift, N, N_shift):

   samples_T = samples.T

   extra_1 = np.full((samples.shape[0],samples.shape[1]), T_shift) 
   extra_2 = np.full((samples_T.shape[1], N.shape[0]), N_shift)

   extra_T = extra_1.T
   extra_N = extra_2.T

   samples_on_initial_polytope = N.dot(T.dot(samples_T) + extra_T) + extra_N 

   return samples_on_initial_polytope


def gmscale(A, iprint, scltol):

   #--------------------------------------------------------------------------------------
   #
   # This function is a translation of a matlab cobra script you may find here:
   # https://github.com/opencobra/cobratoolbox/blob/master/src/analysis/subspaces/gmscale.m 
   #
   # USAGE:
   # cscale, rscale = gmscale(A, iprint, scltol)

   # INPUTS:
   #     A:          is a (m*n) sparse matrix 
   #     scltol:      should be in the range (0.0, 1.0).
   #              Typically `scltol` = 0.9.  A bigger value like 0.99 asks
   #              gmscale` to work a little harder (more passes).   

   # OUTPUTS:
   #    cscale, rscale:    column vectors of column and row scales such that
   #                       `R` (inverse) `A` `C` (inverse) should have entries near 1.0,
   #                        where `R= diag(rscale)`, `C = diag(cscale)`.
   #
   #--------------------------------------------------------------------------------------
   
   m = A.shape[0] ; n = A.shape[1]
   A = np.abs(A)
   A_t = A.T     ##  Attention!!! We will work with the transpose matrix as numpy works row based while matlab is column based
   maxpass = 10
   aratio = 1e+50
   damp = 1e-4
   small = 1e-8
   rscale = np.ones((m,1))
   cscale = np.ones((n,1))

   #----------------------------------------------------
   # Main loop
   #----------------------------------------------------

   for npass in range(maxpass):
      
      rscale[rscale == 0] = 1
      sparse = sp.csr_matrix(1./rscale)
      r_diagonal_elements = np.asarray(sparse.todense())      
      Rinv = diags(r_diagonal_elements[0].T, shape = (m,m)).toarray()
      
      SA = np.dot(Rinv, A)
      SA_T = SA.T
      
      J, I = SA_T.nonzero()
      V = SA.T[SA.T > 0]
      invSA = sp.csr_matrix( (1./V, (J, I)), shape = (n, m) ).T
      
      cmax = np.max(SA, axis = 0)
      cmin = np.max(invSA, axis = 0).data
      cmin = 1./ (cmin + 2.2204e-16)
      
      sratio = np.max(cmax/cmin)
      
      if npass > 0:
         c_product = np.multiply(np.max((np.array((cmin.data, damp*cmax))), axis =0), cmax)
         cscale = np.sqrt(c_product)

      check = aratio*scltol     
      
      if npass >= 2 and sratio >= check:
         break

      if npass == maxpass:
         break
      
      aratio = sratio
      
      # Set new row scales for the next pass.
      
      cscale[cscale == 0] = 1
      sparse = sp.csr_matrix(1./cscale)
      c_diagonal_elements = np.asarray(sparse.todense())
      
      Cinv = diags(c_diagonal_elements[0].T, shape = (n,n)).toarray()      
      SA = np.dot(A, Cinv)
      SA_T = SA.T
      
      J, I = SA_T.nonzero()
      V = SA.T[SA.T > 0]
      invSA = sp.csr_matrix( (1./V, (J, I)), shape = (n, m) ).T
      
      rmax = np.max(SA, axis = 1)
      rmin = np.max(invSA, axis = 1).data
      tmp = rmin + 2.2204e-16
      rmin = 1.0/tmp
      
      r_product = np.multiply(np.max((np.array((rmin.data, damp*rmax))), axis =0), rmax)
      rscale = np.sqrt(r_product)      
      
   #----------------------------------------------------
   # End of main loop
   #----------------------------------------------------

   # Reset column scales so the biggest element in each scaled column will be 1.
   # Again, allow for empty rows and columns.
   rscale[rscale == 0] = 1
   sparse = sp.csr_matrix(1./rscale)
   r_diagonal_elements = np.asarray(sparse.todense())      
   Rinv = diags(r_diagonal_elements[0].T, shape = (m,m)).toarray()

   SA = np.dot(Rinv, A)
   SA_T = SA.T  
   J, I = SA_T.nonzero()
   V = SA.T[SA.T > 0]
   
   cscale = np.max(SA, axis = 0)
   cscale[cscale == 0] = 1
   
   if iprint > 0 :
   
      rmin = np.amin(rscale) ; imin = np.where(rscale == np.amin(rscale))
      cmin = np.amin(cscale) ; jmin = np.where(cscale == np.amin(cscale))
      
      rmax = np.amax(rscale) ; imax = np.where(rscale == np.amax(rscale))
      cmax = np.amax(cscale) ; jmax = np.where(cscale == np.amax(cscale))
      
   return cscale, rscale

def scaled_polytope(full_dim_polytope, cs, rs):
   
   m = rs.shape[0] ; n = cs.shape[0]
   r_diagonal_matrix = diags(1/rs, shape = (m,m)).toarray()
   c_diagonal_matrix = diags(1/cs, shape = (n,n)).toarray()   
   
   P_A = np.dot(r_diagonal_matrix, np.dot(full_dim_polytope.A, c_diagonal_matrix))
   P_b = np.dot(r_diagonal_matrix, full_dim_polytope.b)
   
   return P_A, P_b, c_diagonal_matrix


################################################################################
#                  Classes for the volesti C++ code                            #
################################################################################

# Get classes from the bindings.h file
cdef extern from "bindings.h":

   # The HPolytopeCPP class along with its functions
   cdef cppclass HPolytopeCPP:

      # Initialization
      HPolytopeCPP() except +
      HPolytopeCPP(double *A, double *b, int n_hyperplanes, int n_variables) except +

      # Compute volume
      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed);

      # Random sampling
      double generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, bool boundary, \
         bool cdhr, bool rdhr, bool gaussian, bool set_L, bool accelerated_billiard, bool billiard, bool ball_walk, \
         double a, double L, bool max_ball, double* inner_point, double radius, double* samples);

      # Rounding H-Polytope
      void rounding(char* rounding_method, double* new_A, double* new_b, double* T_matrix, double* shift, double &round_value, \
         bool max_ball, double* inner_point, double radius);
      
      # Rounding svd step 
      double rounding_svd_step(double* new_A, double* new_b, double* T_matrix, double* shift, double* inner_point, double radius)

   # The lowDimPolytopeCPP class along with its functions
   cdef cppclass lowDimHPolytopeCPP:

      # Initialization
      lowDimHPolytopeCPP() except +
      lowDimHPolytopeCPP(double *A, double *b, double *Aeq, double *beq, int n_rows_of_A, int n_cols_of_A, int n_row_of_Aeq, int n_cols_of_Aeq) except +

      # Get full dimensional polytope
      int full_dimensiolal_polytope(double* N_extra_trans, double* shift, double* A_full_extra_trans, double* b_full)

# Lists with the methods supported by volesti for volume approximation and random walk
volume_methods = ["sequence_of_balls".encode("UTF-8"), "cooling_gaussian".encode("UTF-8"), "cooling_balls".encode("UTF-8")]
walk_methods = ["uniform_ball".encode("UTF-8"), "CDHR".encode("UTF-8"), "RDHR".encode("UTF-8"), "gaussian_ball".encode("UTF-8"), \
                "gaussian_CDHR".encode("UTF-8"), "gaussian_RDHR".encode("UTF-8"), "uniform_ball".encode("UTF-8"), "billiard".encode("UTF-8")]
rounding_methods = ["min_ellipsoid".encode("UTF-8"), "svd".encode("UTF-8"), "max_ellipsoid".encode("UTF-8")]

# Build the HPolytope class
cdef class HPolytope:

   cdef HPolytopeCPP polytope_cpp
   cdef double[:,::1] _A
   cdef double[::1] _b

# Set the specs of the class
   def __cinit__(self, double[:,::1] A, double[::1] b):
      self._A = A
      self._b = b
      n_hyperplanes, n_variables = A.shape[0], A.shape[1]
      self.polytope_cpp = HPolytopeCPP(&A[0,0], &b[0], n_hyperplanes, n_variables)

#  This is where the volesti functions are getting their python interface; first the compute_volume() function
   def compute_volume(self, walk_len = 2, epsilon = 0.05, vol_method = "sequence_of_balls", walk_method = "uniform_ball", \
      np.npy_int32 seed=get_time_seed()):

      vol_method = vol_method.encode("UTF-8")
      walk_method = walk_method.encode("UTF-8")

      if vol_method in volume_methods:
         if walk_method in walk_methods:
            return self.polytope_cpp.compute_volume(vol_method, walk_method, walk_len, epsilon, seed)
         else:
            raise Exception('"{}" is not implemented to walk methods. Available methods are: {}'.format(walk_method, walk_methods))
      else:
         raise Exception('"{}" is not implemented to compute volume. Available methods are: {}'.format(vol_method, volume_methods))

# Likewise, the generate_samples() function
   def generate_samples(self, walk_len = 1, number_of_points = 1000, number_of_points_to_burn = 0, boundary = False, cdhr=False, \
      rdhr = False, gaussian = False, set_L = False, accelerated_billiard = True, billiard = False, ball_walk = False, a = 0, \
      radius = 0, inner_point = [], L = 0):

      n_variables = self._A.shape[1]
      cdef double[:,::1] samples = np.zeros((number_of_points,  n_variables), dtype = np.float64, order = "C")
      cdef double[::1] inner_point_for_c = np.asarray(inner_point)
      
      # Check whether the user asks for a certai value of radius; this is of higher priority than having a radius from the corresponding function
      if radius <= 0:        
        max_ball = False
      else:
         max_ball = True
            
      if L <= 0:
         set_L = False
      else:
         set_L = True
      
      self.polytope_cpp.generate_samples(walk_len, number_of_points, number_of_points_to_burn, boundary, cdhr, rdhr, gaussian, set_L, \
                                 accelerated_billiard, billiard, ball_walk, a, L, max_ball, &inner_point_for_c[0], radius, &samples[0,0])
      return np.asarray(samples)      # we need to build a Python function for getting a starting point depending on the polytope

# Special rounding function for the SVD case 
   def rounding_svd(self, scale = None):

      # Get the dimensions of the items about to build
      n_hyperplanes, n_variables = self._A.shape[0], self._A.shape[1]

      # Set the variables of those items; notice that they are all cdef type except of the last one which is about to be used
      # both as a C++ and a Python variable
      cdef double[:,::1] new_A = np.zeros((n_hyperplanes, n_variables), dtype=np.float64, order="C")
      cdef double[::1] new_b = np.zeros(n_hyperplanes, dtype=np.float64, order="C")
      cdef double[:,::1] T_matrix = np.zeros((n_variables, n_variables), dtype=np.float64, order="C")
      cdef double[::1] shift = np.zeros((n_variables), dtype=np.float64, order="C")

      # Get max ball for the initial polytope
      temp_c, radius = get_max_ball(self._A, self._b)
      cdef double[::1] inner_point_for_c = np.asarray(temp_c)
      
      # Build a while loop until for the rounding to converge
      counterrr = 0
      while True:

         counterrr += 1
         check = self.polytope_cpp.rounding_svd_step(&new_A[0,0], &new_b[0], &T_matrix[0,0], &shift[0], &inner_point_for_c[0], radius)
         
         if check < 2.0 and check > 1.0:
            break

         new_temp_c, radius = get_max_ball(new_A, new_b)
         inner_point_for_c = np.asarray(new_temp_c)

      np.save('A_rounded.npy', new_A) ; np.save('b_rounded.npy', new_b)
      np.save('T_rounded.npy', T_matrix) ; np.save('shift_rounded.npy', shift)
      
      try:
         if scale.size != 0:
            T_matrix = np.dot(scale, T_matrix)
      except:
         pass
         
      return np.asarray(new_A), np.asarray(new_b), np.asarray(T_matrix), np.asarray(shift)

     

   
# The rounding() function; like the compute_volume; there are more than one methods for this step
   def rounding(self, rounding_method = 'max_ellipsoid', inner_point = [], radius = 0, scale = None):

      # Get the dimensions of the items about to build
      n_hyperplanes, n_variables = self._A.shape[0], self._A.shape[1]

      # Set the variables of those items; notice that they are all cdef type except of the last one which is about to be used
      # both as a C++ and a Python variable
      cdef double[:,::1] new_A = np.zeros((n_hyperplanes, n_variables), dtype=np.float64, order="C")
      cdef double[::1] new_b = np.zeros(n_hyperplanes, dtype=np.float64, order="C")
      cdef double[:,::1] T_matrix = np.zeros((n_variables, n_variables), dtype=np.float64, order="C")
      cdef double[::1] shift = np.zeros((n_variables), dtype=np.float64, order="C")
      cdef double round_value
      
      cdef double[::1] inner_point_for_c = np.asarray(inner_point)
      
      # Transform the rounding_method variable to UTF-8 coding
      rounding_method = rounding_method.encode("UTF-8")

      # Check whether a max ball has been given
      if radius > 0:
         max_ball = True
      else:
         max_ball = False
      
      # Check whether the rounding method the user asked for, is actually among those volestipy supports
      if rounding_method in rounding_methods:

         self.polytope_cpp.rounding(rounding_method, &new_A[0,0], &new_b[0], &T_matrix[0,0], &shift[0], round_value, max_ball, &inner_point_for_c[0], radius)

         np.save('A_rounded.npy', new_A) ; np.save('b_rounded.npy', new_b)
         np.save('T_rounded.npy', T_matrix) ; np.save('shift_rounded.npy', shift)
         np.save('round_value.npy', np.asarray(round_value))

         try:
            if scale.size != 0:
               T_matrix = np.dot(scale, T_matrix)
         except:
            pass

         return np.asarray(new_A),np.asarray(new_b),np.asarray(T_matrix),np.asarray(shift),np.asarray(round_value)

      else:

         raise Exception('"{}" is not implemented to walk types. Available methods are: {}'.format(rounding_method, rounding_methods))

   @property
   def A(self):
      return np.asarray(self._A)
   @property
   def b(self):
      return np.asarray(self._b)
   @property
   def dimensions(self):
      return self._A.shape[1]

# Build the low_dim_polytope_cpp class
cdef class low_dim_HPolytope:

   cdef lowDimHPolytopeCPP low_dim_polytope_cpp
   cdef double[:,::1] _A
   cdef double[::1] _b
   cdef double [:,::1] _Aeq
   cdef double[::1] _beq

# Set the specs of the class
   def __cinit__(self, double[:,::1] A, double[::1] b, double[:,::1] Aeq, double[::1] beq):

      self._A = A
      self._b = b
      self._Aeq = Aeq
      self._beq = beq
      n_rows_of_A, n_cols_of_A = A.shape[0], A.shape[1]
      n_row_of_Aeq, n_cols_of_Aeq = Aeq.shape[0], Aeq.shape[1]

      # If statements to check whether the user's input is valid for the low_dim_HPolytope class to run
      if n_rows_of_A == b.shape[0]:

         if n_row_of_Aeq == beq.shape[0]:

            if n_cols_of_A == n_cols_of_Aeq:

               # Run the constructor
               self.low_dim_polytope_cpp = lowDimHPolytopeCPP(&A[0,0], &b[0], &Aeq[0,0], &beq[0], n_rows_of_A, n_cols_of_A, n_row_of_Aeq, n_cols_of_Aeq)

            else:
               raise Exception('The number of columns of A equals to "{}" while those of Aeq {}. \
                               A and Aeq need to have the same number of columns'.format(n_cols_of_A, n_cols_of_Aeq))
         else:
            raise Exception('The number of rows of Aeq equals to "{}" while the elements of the beq vector are {}. \
                            The beq vector needs to have length equal to the number of rows of Aeq.'.format(n_row_of_Aeq, beq.shape[0]))
      else:
         raise Exception('The number of rows of A equals to "{}" while the elements of b are {}. \
                         The b vector needs to have length equal to the number of rows of A.'.format(n_rows_of_A, b.shape[0]))

   # The get_full_dimensional_polytope() function(); that needs to run in case the user does not provide volestipy with a full dimensional polytope
   def full_dimensiolal_polytope(self):

      # Get dimensions of the initial S (Aeq) matrix
      m = self._Aeq.shape[0]
      n = self._Aeq.shape[1]
      k = self._A.shape[0]

      # Set the output variables
      # The number of lines in the transpose N (columns in the actual matrix) are at least n-m; but we do not know their exact number
      # So we initialize it with the maximum possible number of lines (n). the same is for the full A transpose matrix
      # Later, we will have to keep their actual dimension and remove these variables with the extra lines
      cdef double[:,::1] N_extra_trans = np.zeros((n, n), dtype=np.float64, order="C")
      cdef double[::1] shift = np.zeros((n), dtype=np.float64, order="C")
      cdef double[:,::1] A_full_extra_trans = np.zeros((n,k), dtype=np.float64, order="C")
      cdef double[::1] b_full = np.zeros((k), dtype=np.float64, order="C")

      # We need to keep the final number of columns of the N / full_A matrices
      cpdef int n_of_cols_in_N

      # Call the C++ class to get the full_dimensional polytope
      n_of_cols_in_N = self.low_dim_polytope_cpp.full_dimensiolal_polytope(&N_extra_trans[0,0], &shift[0], &A_full_extra_trans[0,0], &b_full[0])

      # Get a matrix with exactly the number of lines and columns that N expands to and delete the one with the extra columns
      N = np.zeros((n, n_of_cols_in_N), dtype=np.float64, order="C")
      for i in range(n):
         for j in range(n_of_cols_in_N):
            N[i,j] = np.asarray(N_extra_trans[j,i])
      del N_extra_trans

      # Likewise, for the A matrix of the full dimensional polytope
      A_full = np.zeros((k, n_of_cols_in_N), dtype=np.float64, order="C")
      for i in range(k):
         for j in range(n_of_cols_in_N):
            A_full[i,j] = np.asarray(A_full_extra_trans[j,i])
      del A_full_extra_trans

      # Finally, we need to build an HP object for the full dumensional polytope we got
      full_dimensional_polytope = HPolytope(A_full,b_full)

      # Print all the output of the function in .npy files
      np.save('A_full_dim.npy', A_full) ; np.save('b_full_dim.npy', b_full)
      np.save('N_full_dim.npy', N) ; np.save('shift_full_dim.npy',shift)

      # Delete all non-needed vars
      del A_full
      del b_full

      # Return a tuple whith the full dimensional HPolytope object in the first position ([0]) the N matrix and the shift vector
      return full_dimensional_polytope, np.asarray(N), np.asarray(shift)


   @property
   def A(self):
      return np.asarray(self._A)
   @property
   def b(self):
      return np.asarray(self._b)
   @property
   def Aeq(self):
      return np.asarray(self._Aeq)
   @property
   def beq(self):
      return np.asarray(self._beq)
   @property
   def dimensions(self):
      return self._A.shape[1]
   
