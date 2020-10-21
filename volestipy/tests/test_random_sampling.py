#!/usr/bin/python3.6

import numpy as np
from volestipy import HPolytope


if __name__ == "__main__":

    dim = 3
    A = np.zeros((2*dim, dim), dtype=np.float)
    A[0:dim] = np.eye(dim)
    A[dim:] -=  np.eye(dim,dim, dtype=np.float)
    b = np.ones(2*dim, dtype=np.float)
    p = HPolytope(A,b)

    print("*** This is a test for the random_sampling() function of the volestipy library ***\n")
    print("The polytope:")
    print("Dimensions: {}".format(p.dimensions))
    print("A matrix:")
    print(p.A)
    print("b vector:")
    print(p.b)

    samples = p.generate_samples(walk_len = 5, number_of_points = 80000, number_of_points_to_burn = 50, boundary = False, cdhr = True, rdhr = False, 
              gaussian = False, set_L = False, billiard = False, ball_walk = False, a = 0, L = 0)


    print("\n >> This is the output for the CDHR random sampling algorithm using the boundary option <<\n")
    print("Samples:")
    print(samples) 


