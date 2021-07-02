// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>
#include "misc.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"

#include "random_walks/random_walks.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"

#include "preprocess/max_inscribed_ellipsoid_rounding.hpp"

#include "known_polytope_generators.h"

template <typename NT>
NT factorial(NT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename NT>
void test_values(NT volume, NT expected, NT exact)
{
    std::cout << "Computed volume " << volume << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
    std::cout << "Relative error (expected) = "
              << std::abs((volume-expected)/expected) << std::endl;
    std::cout << "Relative error (exact) = "
              << std::abs((volume-exact)/exact) << std::endl;
    CHECK((std::abs((volume - exact)/exact) < 0.2 || 
            std::abs((volume - expected)/expected) < 0.00001));
}

template <class Polytope>
void rounding_test(Polytope &HP,
                 double const& expectedBall,
                 double const& expectedCDHR,
                 double const& expectedRDHR,
                 double const& expectedBilliard,
                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    int d = HP.dimension();

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;
    RNGType rng(d);

    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();
    std::tuple<MT, VT, NT> res = max_inscribed_ellipsoid_rounding<MT, VT, NT>(HP, InnerBall.first);

    // Setup the parameters
    int walk_len = 10;
    NT e = 0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;


    //TODO: low accuracy in high dimensions
    // NT volume = res.second * volume_cooling_balls<BallWalk, RNGType>(HP, e, walk_len);
    // test_values(volume, expectedBall, exact);

    NT volume = std::get<2>(res) * volume_cooling_balls<CDHRWalk, RNGType>(HP, e, walk_len).second;
    test_values(volume, expectedCDHR, exact);

    volume = std::get<2>(res) * volume_cooling_balls<RDHRWalk, RNGType>(HP, e, 2*walk_len).second;
    test_values(volume, expectedRDHR, exact);

    volume = std::get<2>(res) * volume_cooling_balls<BilliardWalk, RNGType>(HP, e, walk_len).second;
    test_values(volume, expectedBilliard, exact);
}




template <typename NT>
void call_test_skinny_cubes() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing rounding of H-skinny_cube5" << std::endl;
    P = generate_skinny_cube<Hpolytope>(5);
    rounding_test(P, 0, 3070.64, 3188.25, 3140.6, 3200.0);

    std::cout << "\n--- Testing rounding of H-skinny_cube10" << std::endl;

    P = generate_skinny_cube<Hpolytope>(10);
    rounding_test(P, 0, 122550, 108426, 105003.0, 102400.0);

    std::cout << "\n--- Testing rounding of H-skinny_cube20" << std::endl;
    P = generate_skinny_cube<Hpolytope>(20);
    rounding_test(P, 0,
                  8.26497 * std::pow(10,7),
                  8.94948e+07,
                  1.09218e+08,
                  104857600.0);
}


TEST_CASE("round_skinny_cube") {
    call_test_skinny_cubes<double>();
}
