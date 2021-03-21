// VolEsti (volume computation and sampling library)

// Copyright (c) 2020 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_RANDOM_WALKS_HPP
#define RANDOM_WALKS_RANDOM_WALKS_HPP

#include "random_walks/boundary_cdhr_walk.hpp"
#include "random_walks/boundary_rdhr_walk.hpp"
#include "random_walks/gaussian_ball_walk.hpp"
#include "random_walks/gaussian_cdhr_walk.hpp"
#include "random_walks/gaussian_rdhr_walk.hpp"
#include "random_walks/uniform_ball_walk.hpp"
#include "random_walks/uniform_billiard_walk.hpp"
#include "random_walks/uniform_cdhr_walk.hpp"
#include "random_walks/uniform_rdhr_walk.hpp"
#include "random_walks/uniform_dikin_walk.hpp"
#include "random_walks/uniform_john_walk.hpp"
#include "random_walks/uniform_vaidya_walk.hpp"
#include "random_walks/uniform_accelerated_billiard_walk.hpp"
#include "random_walks/exponential_hamiltonian_monte_carlo_walk.hpp"
#include "random_walks/exponential_hmc_leapfrog.hpp"
#ifndef VOLESTIPY
    #include "random_walks/hamiltonian_monte_carlo_walk.hpp"
    #include "random_walks/langevin_walk.hpp"
#endif

#endif // RANDOM_WALKS_RANDOM_WALKS_HPP
