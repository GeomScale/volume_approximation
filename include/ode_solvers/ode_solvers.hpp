// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

/*
  Each solver solves an ODE of the form d^k x/dt^k = F(x, t)
  The vector x can be written as x = (x_1, ...., x_n) where each
  x_i is conditioned on a convex polytope K_i with bounary reflections
  on the boundary of K_i for each 1 <= i <= n. Each sub-state is determined
  by the oracle function F_i(x, t) which has range over x = (x_1, ..., x_n).

  Some general parameter notations for the solvers

  Templates
      1. Polytope: The polytope type (H-Polytope, V-Polytope, Z-Polytope) K_i
      2. NT: Number type (double, float)
      3. Point: Point type
      4. func: Function type for the oracle collection K_i
      5. bfunc: Basis function type (e.g. polynomial) for non-linear trajectory
         methods (such as collocation)

  Variables
      1. Ks: A vector of domains K_i for 1 <= i <= n
      2. xs: A vector of substates x_i for 1 <= i <= n
      3. xs_prev: The previous state of the solver
      4. F: A functor of oracles F_i for 1 <= i <= n
      5. eta: Step size
      6. t: Temporal variable

  TODO:
      1. Change datastructure of boundaries (std::vector<Polytope*>)

*/

#include <iostream>
#include <cmath>
#include <functional>
#include <vector>

#include "ode_solvers/euler.hpp"
#include "ode_solvers/runge_kutta.hpp"
#include "ode_solvers/leapfrog.hpp"
#include "ode_solvers/richardson_extrapolation.hpp"
#include "ode_solvers/oracle_functors.hpp"
#include "ode_solvers/randomized_midpoint.hpp"
#include "ode_solvers/generalized_leapfrog.hpp"

#ifndef DISABLE_NLP_ORACLES
#include "ode_solvers/collocation.hpp"
#include "ode_solvers/basis.hpp"
#include "ode_solvers/integral_collocation.hpp"
#endif

#ifndef ODE_SOLVERS_ODE_SOLVERS_HPP
#define ODE_SOLVERS_ODE_SOLVERS_HPP

#endif
