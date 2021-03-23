// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

#include "Eigen/Eigen"

#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "sampling/sampling.hpp"
#include "generators/known_polytope_generators.h"
#include "diagnostics/multivariate_psrf.hpp"


template <typename NT>
void run_main() 
{
    typedef Cartesian<NT>    Kernel;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef typename Kernel::Point    Point;
    typedef HPolytope <Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    int dim = 50;

    Hpolytope P = generate_cube<Hpolytope>(dim, false);
    std::list<Point> randPoints;
    RNGType rng(dim);
    unsigned int walkL = 1, numpoints = 1000, nburns = 0;
    NT variance = 1.0;
    NT a = NT(1) / (NT(2) * variance);

    Point c(dim), StartingPoint(dim);
    P.set_InnerBall(std::pair<Point,NT>(Point(dim), 1.0));
    c = GetDirection<Point>::apply(dim, rng, false);
    gaussian_sampling<GaussianHamiltonianMonteCarloExactWalk>(randPoints, P, rng, walkL, numpoints, a,
                                                              StartingPoint, nburns);
    MT samples;
    samples.resize(dim, numpoints);
    int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++) 
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }

    std::cerr << "PSRF: " <<  multivariate_psrf<NT, VT, MT>(samples) << std::endl;
}

int main() {
  run_main<double>();
  return 0;
}
