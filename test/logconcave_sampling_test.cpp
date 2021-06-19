// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <unistd.h>
#include <string>
#include <typeinfo>
#include <fstream>
#include <thread>
#include <atomic>
#include <tuple>

#include "doctest.h"
#include "Eigen/Eigen"

#include "ode_solvers.hpp"
#include "diagnostics/diagnostics.hpp"

#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"

#include "generators/known_polytope_generators.h"
#include "generators/h_polytopes_generator.h"
#include "generators/convex_bodies_generator.h"

#include "preprocess/svd_rounding.hpp"
#include "misc/misc.h"

#include "matrix_operations/EigenvaluesProblems.h"
#include "SDPAFormatManager.h"
#include "convex_bodies/spectrahedra_new/newSpectrahedron.h"

template <typename NT>
struct SimulationStats {
    std::string method;
    unsigned int walk_length;
    unsigned int min_ess = 0;
    NT max_psrf = NT(0);
    NT time_per_draw = NT(0);
    NT time_per_independent_sample = NT(0);
    NT average_acceptance_log_prob = NT(0);
    NT average_number_of_reflections = NT(0);
    NT step_size = NT(0);


    friend std::ostream& operator<< (std::ostream& out, const SimulationStats &stats) {
        out << stats.method << ","
            << stats.walk_length << ","
            << stats.min_ess << ","
            << stats.max_psrf << ","
            << stats.time_per_draw << ","
            << stats.time_per_independent_sample << ","
            << stats.average_acceptance_log_prob << ","
            << stats.average_number_of_reflections << ","
            << stats.step_size
            << std::endl;
        return out;
    }
};

struct InnerBallFunctor {

  // Gaussian density centered at the inner ball center
  template <
      typename NT,
      typename Point
  >
  struct parameters {
    unsigned int order;
    NT L; // Lipschitz constant for gradient
    NT m; // Strong convexity constant
    NT kappa; // Condition number
    NT R0;
    NT sigma;
    Point x0;

    parameters(Point x0_, NT R0_) : order(2), L(1), m(1), kappa(1), x0(x0_), R0(R0_), sigma(1) {};

  };

  template
  <
      typename Point
  >
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT, Point> &params;

    GradientFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
      if (i == params.order - 1) {
        Point y = (-1.0 / pow(params.sigma, 2)) * (xs[0] - params.x0);
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }

  };

  template
  <
    typename Point
  >
  struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT, Point> &params;

    FunctionFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    NT operator() (Point const& x) const {
      Point y = x - params.x0;
      return 1.0 / (2 * pow(params.sigma, 2)) * y.dot(y);
    }

  };

};

struct CustomFunctor {

  // Custom density with neg log prob equal to || x ||^2 + 1^T x
  template <
      typename NT
  >
  struct parameters {
    unsigned int order;
    NT L; // Lipschitz constant for gradient
    NT m; // Strong convexity constant
    NT kappa; // Condition number

    parameters() : order(2), L(2), m(2), kappa(1) {};

    parameters(unsigned int order_) :
      order(order),
      L(2),
      m(2),
      kappa(1)
    {}
  };

  template
  <
      typename Point
  >
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT> params;

    GradientFunctor() {};

    // The index i represents the state vector index
    Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
      if (i == params.order - 1) {
        Point y = (-1.0) * Point::all_ones(xs[0].dimension());
        y = y + (-2.0) * xs[0];
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }

  };

  template
  <
    typename Point
  >
  struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT> params;

    FunctionFunctor() {};

    // The index i represents the state vector index
    NT operator() (Point const& x) const {
      return x.dot(x) + x.sum();
    }

  };

};

template <typename NT, typename VT, typename MT>
NT check_interval_psrf(MT &samples, NT target=NT(1.2)) {
    NT max_psrf = NT(0);
    VT intv_psrf = interval_psrf<VT, NT, MT>(samples);
    unsigned int d = intv_psrf.rows();
    for (unsigned int i = 0; i < d; i++) {
        CHECK(intv_psrf(i) < target);
        if (intv_psrf(i) > max_psrf) max_psrf = intv_psrf(i);
    }
    return max_psrf;
}

template <typename Sampler, typename RandomNumberGenerator, typename NT, typename Point>
void check_ergodic_mean_norm(
    Sampler &sampler,
    RandomNumberGenerator &rng,
    Point &mean,
    unsigned int dim,
    int n_samples=1500,
    int skip_samples=750,
    NT target=NT(0),
    NT tol=5e-1) {

  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < n_samples; i++) {
    sampler.apply(rng, 1);
    if (i >= skip_samples) {
      mean = mean + sampler.x;
    }

    #ifdef VOLESTI_DEBUG
      std::cout << sampler.x.getCoefficients().transpose() << std::endl;
    #endif
  }

  auto stop = std::chrono::high_resolution_clock::now();

  long ETA = (long) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

  mean = (1.0 / (n_samples - skip_samples)) * mean;

  NT error = abs(NT(mean.dot(mean)) - target);

  if (target != NT(0)) error /= abs(target);

  std::cout << "Dimensionality: " << dim << std::endl;
  std::cout << "Target ergodic mean norm: " << target << std::endl;
  std::cout << "Error (relative if possible) after " << n_samples << " samples: " << error << std::endl;
  std::cout << "ETA (us): " << ETA << std::endl << std::endl;

  CHECK(error < tol);

}

template <typename NT>
void benchmark_hmc(bool truncated) {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef HPolytope<Point> Hpolytope;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
    typedef CustomFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef CustomFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef LeapfrogODESolver<Point, NT, Hpolytope, NegativeGradientFunctor> Solver;

    NegativeGradientFunctor F;
    NegativeLogprobFunctor f;
    RandomNumberGenerator rng(1);
    unsigned int dim_min = 1;
    unsigned int dim_max = 100;
    int n_samples = 1000;
    long ETA;
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;

    for (unsigned int dim = dim_min; dim <= dim_max; dim++) {
      Point x0(dim);
      HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);
      if (truncated) {
        Hpolytope P = generate_cube<Hpolytope>(dim, false);
        HamiltonianMonteCarloWalk::Walk<Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
        hmc(&P, x0, F, f, hmc_params);
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n_samples; i++) hmc.apply(rng, 1);
        stop = std::chrono::high_resolution_clock::now();

      } else {

        HamiltonianMonteCarloWalk::Walk<Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
        hmc(NULL, x0, F, f, hmc_params);

        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n_samples; i++) hmc.apply(rng, 1);
        stop = std::chrono::high_resolution_clock::now();
      }

      ETA = (long) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
      std::cout << ETA << std::endl;
    }

}

template <typename NT>
void test_hmc() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef HPolytope<Point> Hpolytope;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
    typedef IsotropicQuadraticFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef IsotropicQuadraticFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef LeapfrogODESolver<Point, NT, Hpolytope, NegativeGradientFunctor> Solver;

    IsotropicQuadraticFunctor::parameters<NT> params;
    params.order = 2;

    NegativeGradientFunctor F(params);
    NegativeLogprobFunctor f(params);

    RandomNumberGenerator rng(1);
    unsigned int dim = 10;
    HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);
    Hpolytope P = generate_cube<Hpolytope>(dim, false);
    Point x0(dim);

    HamiltonianMonteCarloWalk::Walk
      <Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      hmc(&P, x0, F, f, hmc_params);

    Point mean(dim);
    check_ergodic_mean_norm(hmc, rng, mean, dim, 75000, 37500, NT(0));
}


template <typename NT>
void test_uld() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef HPolytope<Point> Hpolytope;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
    typedef IsotropicQuadraticFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef IsotropicQuadraticFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef RandomizedMipointSDESolver<Point, NT, Hpolytope, NegativeGradientFunctor, RandomNumberGenerator> Solver;

    IsotropicQuadraticFunctor::parameters<NT> params;
    params.order = 2;
    params.alpha = NT(1);

    NegativeGradientFunctor F(params);
    NegativeLogprobFunctor f(params);

    RandomNumberGenerator rng(1);
    unsigned int dim = 5;
    UnderdampedLangevinWalk::parameters<NT, NegativeGradientFunctor> uld_params(F, dim);
    Hpolytope P = generate_cube<Hpolytope>(dim, false);
    Point x0(dim);

    UnderdampedLangevinWalk::Walk
      <Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      uld(&P, x0, F, f, uld_params);

    Point mean(dim);
    check_ergodic_mean_norm(uld, rng, mean, dim, 75000, 37500, NT(0));

}

template <typename NT, typename Polytope>
std::vector<SimulationStats<NT>> benchmark_polytope_sampling(
    Polytope &P,
    NT eta=NT(-1),
    unsigned int walk_length=3,
    bool rounding=false,
    bool centered=false,
    bool warmstart=true,
    unsigned int max_draws=80000,
    unsigned int num_burns=20000) {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef boost::mt19937 RNGType;
    typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
    typedef InnerBallFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef InnerBallFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef LeapfrogODESolver<Point, NT, Polytope, NegativeGradientFunctor> Solver;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    SimulationStats<NT> rdhr_stats;
    SimulationStats<NT> hmc_stats;

    std::pair<Point, NT> inner_ball;
    if (centered) {
        inner_ball.first = Point(P.dimension());
        inner_ball.second = NT(1); // dummy radius (not correct one)
    } else {
        inner_ball = P.ComputeInnerBall();
    }

    // Random number generator
    RandomNumberGenerator rng(1);

    // Chebyshev center
    Point x0 = inner_ball.first;
    NT R0 = inner_ball.second;
    unsigned int dim = x0.dimension();

    if (rounding) {
        std::cout << "SVD Rounding" << std::endl;
        svd_rounding<AcceleratedBilliardWalk, MT, VT>(P, inner_ball, walk_length, rng);
    }

    // Declare oracles
    InnerBallFunctor::parameters<NT, Point> params(x0, R0);

    NegativeGradientFunctor F(params);
    NegativeLogprobFunctor f(params);

    GaussianRDHRWalk::Walk<Polytope, RandomNumberGenerator> gaussian_walk(P, x0, params.L, rng);

    int max_actual_draws = max_draws - num_burns;
    unsigned int min_ess = 0;

    MT samples;
    samples.resize(dim, max_actual_draws);
    NT ETA;
    NT max_psrf;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;

    if (warmstart) {
        std::cout << "Gaussian Hit and Run" << std::endl;

        std::cout << "Burn-in" << std::endl;

        for (unsigned int i = 0; i < num_burns; i++) {
          if (i % 1000 == 0) std::cout << ".";
          gaussian_walk.apply(P, x0, params.L, walk_length, rng);
          // std::cout << x0.getCoefficients() << std::endl;
        }

        std::cout << std::endl;
        std::cout << "Sampling" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        for (unsigned int i = 0; i < max_actual_draws; i++) {
          gaussian_walk.apply(P, x0, params.L, walk_length, rng);
          samples.col(i) = x0.getCoefficients();
          if (i % 1000 == 0 && i > 0) std::cout << ".";
        }
        stop = std::chrono::high_resolution_clock::now();

        ETA = (NT) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

        std::cout << std::endl;
        print_diagnostics<NT, VT, MT>(samples, min_ess, std::cout);
        std::cout << "Average time per sample: " << ETA / max_actual_draws << "us" << std::endl;
        std::cout << "Average time per independent sample: " << ETA / min_ess << "us" << std::endl;
        std::cout << std::endl;

        max_psrf = check_interval_psrf<NT, VT, MT>(samples);

        rdhr_stats.method = "RDHR";
        rdhr_stats.walk_length = walk_length;
        rdhr_stats.min_ess = min_ess;
        rdhr_stats.max_psrf = max_psrf;
        rdhr_stats.time_per_draw = ETA / max_actual_draws;
        rdhr_stats.time_per_independent_sample = ETA / min_ess;
    }

    HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);

    HamiltonianMonteCarloWalk::Walk
      <Point, Polytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      hmc(&P, x0, F, f, hmc_params);

    min_ess = 0;

    std::cout << "Hamiltonian Monte Carlo (Gaussian Density)" << std::endl;

    if (eta > 0) hmc.solver->eta = eta;

    std::cout << "Burn-in" << std::endl;

    for (unsigned int i = 0; i < num_burns; i++) {
      if (i % 1000 == 0) std::cout << ".";
      hmc.apply(rng, walk_length);
    }

    hmc.disable_adaptive();
    std::cout << std::endl;
    std::cout << "Sampling" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < max_actual_draws; i++) {
      hmc.apply(rng, walk_length);
      samples.col(i) = hmc.x.getCoefficients();
      if (i % 1000 == 0 && i > 0) std::cout << ".";
    }
    stop = std::chrono::high_resolution_clock::now();

    ETA = (NT) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

    std::cout << std::endl;
    print_diagnostics<NT, VT, MT>(samples, min_ess, std::cout);
    std::cout << "min ess " << min_ess << "us" << std::endl;
    std::cout << "Average time per sample: " << ETA / max_actual_draws << "us" << std::endl;
    std::cout << "Average time per independent sample: " << ETA / min_ess << "us" << std::endl;
    std::cout << "Average number of reflections: " <<
        (1.0 * hmc.solver->num_reflections) / hmc.solver->num_steps << std::endl;
    std::cout << "Step size (final): " << hmc.solver->eta << std::endl;
    std::cout << "Discard Ratio: " << hmc.discard_ratio << std::endl;
    std::cout << "Average Acceptance Probability: " << exp(hmc.average_acceptance_log_prob) << std::endl;
    std::cout << std::endl;

    max_psrf = check_interval_psrf<NT, VT, MT>(samples);

    hmc_stats.method = "HMC";
    hmc_stats.walk_length = walk_length;
    hmc_stats.min_ess = min_ess;
    hmc_stats.max_psrf = max_psrf;
    hmc_stats.time_per_draw = ETA / max_actual_draws;
    hmc_stats.time_per_independent_sample = ETA / min_ess;
    hmc_stats.average_number_of_reflections = (1.0 * hmc.solver->num_reflections) / hmc.solver->num_steps;
    hmc_stats.step_size = hmc.solver->eta;
    hmc_stats.average_acceptance_log_prob  = exp(hmc.average_acceptance_log_prob);

    return std::vector<SimulationStats<NT>>{rdhr_stats, hmc_stats};
}


template <typename NT, typename Polytope>
std::vector<SimulationStats<NT>> benchmark_spectrahedron_sampling(
    Polytope &P,
    NT eta=NT(-1),
    unsigned int walk_length=3,
    bool centered=false,
    unsigned int max_draws=80000,
    unsigned int num_burns=20000) {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef boost::mt19937 RNGType;
    typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
    typedef InnerBallFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef InnerBallFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef LeapfrogODESolver<Point, NT, Polytope, NegativeGradientFunctor> Solver;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    SimulationStats<NT> rdhr_stats;
    SimulationStats<NT> hmc_stats;

    std::pair<Point, NT> inner_ball;
    if (centered) {
        inner_ball.first = Point(P.dimension());
        inner_ball.second = NT(1); // dummy radius (not correct one)
    }

    // Random number generator
    RandomNumberGenerator rng(1);

    // Chebyshev center
    Point x0 = inner_ball.first;
    NT R0 = inner_ball.second;
    unsigned int dim = x0.dimension();
    std::cout<<"dim = "<<dim<<std::endl;


    // Declare oracles
    InnerBallFunctor::parameters<NT, Point> params(x0, R0);

    NegativeGradientFunctor F(params);
    NegativeLogprobFunctor f(params);

    GaussianRDHRWalk::Walk<Polytope, RandomNumberGenerator> gaussian_walk(P, x0, params.L, rng);

    int max_actual_draws = max_draws - num_burns;
    unsigned int min_ess = 0;

    MT samples;
    samples.resize(dim, max_actual_draws);
    NT ETA;
    NT max_psrf;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;

    HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);

    HamiltonianMonteCarloWalk::Walk
      <Point, Polytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      hmc(&P, x0, F, f, hmc_params);

    min_ess = 0;

    std::cout << "Hamiltonian Monte Carlo (Gaussian Density)" << std::endl;

    if (eta > 0) hmc.solver->eta = eta;

    std::cout << "Burn-in" << std::endl;

    for (unsigned int i = 0; i < num_burns; i++) {
      if (i % 1000 == 0) std::cout << ".";
      hmc.apply(rng, walk_length);
    }

    hmc.disable_adaptive();
    std::cout << std::endl;
    std::cout << "Sampling" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < max_actual_draws; i++) {
      hmc.apply(rng, walk_length);
      samples.col(i) = hmc.x.getCoefficients();
      if (i % 1000 == 0 && i > 0) std::cout << ".";
    }
    stop = std::chrono::high_resolution_clock::now();

    ETA = (NT) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

    std::cout << std::endl;
    print_diagnostics<NT, VT, MT>(samples, min_ess, std::cout);
    std::cout << "Average time per sample: " << ETA / max_actual_draws << "us" << std::endl;
    std::cout << "Average time per independent sample: " << ETA / min_ess << "us" << std::endl;
    std::cout << "Average number of reflections: " <<
        (1.0 * hmc.solver->num_reflections) / hmc.solver->num_steps << std::endl;
    std::cout << "Step size (final): " << hmc.solver->eta << std::endl;
    std::cout << "Discard Ratio: " << hmc.discard_ratio << std::endl;
    std::cout << "Average Acceptance Probability: " << exp(hmc.average_acceptance_log_prob) << std::endl;
    std::cout << std::endl;

    max_psrf = check_interval_psrf<NT, VT, MT>(samples);

    hmc_stats.method = "HMC";
    hmc_stats.walk_length = walk_length;
    hmc_stats.min_ess = min_ess;
    hmc_stats.max_psrf = max_psrf;
    hmc_stats.time_per_draw = ETA / max_actual_draws;
    hmc_stats.time_per_independent_sample = ETA / min_ess;
    hmc_stats.average_number_of_reflections = (1.0 * hmc.solver->num_reflections) / hmc.solver->num_steps;
    hmc_stats.step_size = hmc.solver->eta;
    hmc_stats.average_acceptance_log_prob  = exp(hmc.average_acceptance_log_prob);

    return std::vector<SimulationStats<NT>>{rdhr_stats, hmc_stats};
}



template <typename NT, typename Polytope>
void benchmark_multi_chain(
    unsigned int num_chains,
    Polytope &P,
    NT eta=NT(-1),
    unsigned int walk_length=3,
    bool rounding=true,
    unsigned int max_draws=80000,
    unsigned int num_burns=20000) {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef boost::mt19937 RNGType;
    typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    std::pair<Point, NT> inner_ball = P.ComputeInnerBall();

    // Random number generator
    RandomNumberGenerator rng(1);

    // Chebyshev center
    Point x0 = inner_ball.first;
    NT R0 = inner_ball.second;
    unsigned int dim = x0.dimension();

    if (rounding) {
        std::cout << "SVD Rounding" << std::endl;
        svd_rounding<AcceleratedBilliardWalk, MT, VT>(P, inner_ball, walk_length, rng);
    }

    for (unsigned int i = 0; i < num_chains; i++) {
        std::thread chain(benchmark_polytope_sampling<NT, Polytope>, P, eta, walk_length, false, max_draws, num_burns, num_chains);
        chain.join();
    }

}

template <typename NT, typename Polytope, typename Point>
void benchmark_polytope_linear_program_optimization(
    Point &coeffs,
    Polytope &P,
    NT eta=NT(-1),
    unsigned int walk_length=3,
    bool rounding=true,
    unsigned int max_draws=80000,
    unsigned int num_burns=20000) {
    typedef Cartesian<NT>    Kernel;
    typedef std::vector<Point> pts;
    typedef boost::mt19937 RNGType;
    typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
    typedef LinearProgramFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef LinearProgramFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef OptimizationFunctor::GradientFunctor<Point, NegativeLogprobFunctor,
        NegativeGradientFunctor> NegativeGradientOptimizationFunctor;
    typedef OptimizationFunctor::FunctionFunctor<Point, NegativeLogprobFunctor,
        NegativeGradientFunctor> NegativeLogprobOptimizationFunctor;

    typedef LeapfrogODESolver<Point, NT, Polytope,  NegativeGradientOptimizationFunctor> Solver;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    std::pair<Point, NT> inner_ball = P.ComputeInnerBall();

    // Random number generator
    RandomNumberGenerator rng(1);

    // Chebyshev center
    Point x0 = inner_ball.first;
    NT R0 = inner_ball.second;
    unsigned int dim = x0.dimension();

    if (rounding) {
        std::cout << "SVD Rounding" << std::endl;
        svd_rounding<AcceleratedBilliardWalk, MT, VT>(P, inner_ball, walk_length, rng);
    }

    // Declare oracles for LP
    LinearProgramFunctor::parameters<NT, Point> lp_params(coeffs);

    NegativeGradientFunctor F_lp(lp_params);
    NegativeLogprobFunctor f_lp(lp_params);

    // Declare optimization oracles
    OptimizationFunctor::parameters
        <NT, NegativeLogprobFunctor, NegativeGradientFunctor>
        opt_params(1, x0.dimension(), f_lp, F_lp);

    NegativeLogprobOptimizationFunctor f(opt_params);
    NegativeGradientOptimizationFunctor F(opt_params);

    GaussianRDHRWalk::Walk<Polytope, RandomNumberGenerator> gaussian_walk(P, x0, lp_params.L, rng);
    int n_warmstart_samples = 100;

    for (int i = 0; i < n_warmstart_samples; i++) {
        gaussian_walk.apply(P, x0, lp_params.L, walk_length, rng);
    }

    HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientOptimizationFunctor> hmc_params(F, dim);

    HamiltonianMonteCarloWalk::Walk
      <Point, Polytope, RandomNumberGenerator, NegativeGradientOptimizationFunctor, NegativeLogprobOptimizationFunctor, Solver>
      hmc(&P, x0, F, f, hmc_params);

    int max_actual_draws = max_draws - num_burns;
    unsigned int min_ess = 0;

    Point minimum = x0;

    if (eta > 0) hmc.solver->eta = eta;

    std::cout << "Burn-in" << std::endl;

    for (unsigned int i = 0; i < num_burns; i++) {
      if (i % 1000 == 0) std::cout << ".";
      hmc.apply(rng, walk_length);
    }

    std::cout << std::endl;
    std::cout << "Optimizing" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    for (unsigned int j = 0; j < (unsigned int) 4 * ceil(sqrt(dim)); j++) {
        std::cout << "Temperature " << opt_params.T << std::endl;
        for (unsigned int i = 0; i < max_actual_draws; i++) {
            hmc.apply(rng, walk_length);
            if (f_lp(minimum) >= f_lp(hmc.x)) {
                minimum = hmc.x;
                std::cout << "Current value: " << f_lp(minimum) << std::endl;
            }
            if (i % 1000 == 0 && i > 0) std::cout << ".";
        }
        opt_params.update_temperature();
    }
    auto stop = std::chrono::high_resolution_clock::now();

    NT ETA = (NT) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

    std::cout << "LP Value: " << f_lp(minimum) << std::endl;
    std::cout << "Point: " << minimum.getCoefficients().transpose() << std::endl;
}


template <typename Polytope, typename NT>
Polytope read_polytope(std::string filename) {
    std::ifstream inp;
    std::vector<std::vector<NT> > Pin;
    inp.open(filename,std::ifstream::in);
    read_pointset(inp, Pin);
    Polytope P(Pin);
    return P;
}

template <typename NT, typename Point>
std::pair<Point, NT> read_inner_ball_(std::string filename) {
    std::ifstream inp;
    inp.open(filename,std::ifstream::in);
    return read_inner_ball<NT, Point>(inp);
}

template <typename NT>
void call_test_hmc() {
  std::cout << "--- Testing Hamiltonian Monte Carlo" << std::endl;
  test_hmc<NT>();
}

template <typename NT>
void call_test_uld() {
  std::cout << "--- Testing Underdamped Langevin Diffusion" << std::endl;
  test_uld<NT>();
}

template <typename NT>
void call_test_benchmark_hmc(bool truncated) {
  benchmark_hmc<NT>(truncated);
}

template <typename NT>
void call_test_benchmark_polytopes_grid_search() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef boost::mt19937 RNGType;

    std::cout << " --- Grid search on polytopes " << std::endl;

    std::vector<SimulationStats<NT>> results;

    std::vector<std::tuple<Hpolytope, std::string, bool>> polytopes{
       std::make_tuple(generate_skinny_cube<Hpolytope>(100, false), "100_skinny_cube", false),
       std::make_tuple(generate_cross<Hpolytope>(10, false), "10_cross", false),
       std::make_tuple(generate_simplex<Hpolytope>(100, false), "100_simplex", false),
       std::make_tuple(generate_cube<Hpolytope>(100, false), "100_cube", false),
       std::make_tuple(generate_prod_simplex<Hpolytope>(50, false), "50_prod_simplex", false),
       std::make_tuple(generate_birkhoff<Hpolytope>(10), "10_birkhoff", false),
       std::make_tuple(read_polytope<Hpolytope, NT>("./metabolic_full_dim/polytope_iAB_RBC_283.ine"), "iAB_RBC_283", true),
       std::make_tuple(read_polytope<Hpolytope, NT>("./metabolic_full_dim/polytope_iAT_PLT_636.ine"), "iAT_PLT_636", true),
       std::make_tuple(read_polytope<Hpolytope, NT>("./metabolic_full_dim/polytope_e_coli.ine"), "e_coli", true),
       std::make_tuple(read_polytope<Hpolytope, NT>("./metabolic_full_dim/polytope_recon2.ine"), "recon2", true)
    };

    Hpolytope P;
    std::string name;
    std::ofstream outfile;
    NT step_size = 0;
    std::pair<Point, NT> inner_ball;

    for (std::tuple<Hpolytope, std::string, bool> polytope_tuple : polytopes) {
        P = std::get<0>(polytope_tuple);
        name = std::get<1>(polytope_tuple);
        std::cout << name << std::endl;
        outfile.open("results_" + name + "_new.txt");
        P.normalize();
        inner_ball = P.ComputeInnerBall();
        step_size = inner_ball.second / 10;
        for (unsigned int walk_length = 500; walk_length <= P.dimension(); walk_length += P.dimension() / 10) {
            results = benchmark_polytope_sampling<NT, Hpolytope>(P, step_size, walk_length, false, std::get<2>(polytope_tuple), false);
            outfile << results[0];
            outfile << results[1];
        }
        outfile.close();
    }

}


template <typename NT>
void call_test_benchmark_spectrahedra_grid_search() {

    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point    Point;
    typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;
    typedef boost::mt19937 RNGType;

    std::cout << " --- Grid search on spectrahedra " << std::endl;

    std::vector<SimulationStats<NT>> results;
    std::vector<std::tuple<SPECTRAHEDRON, std::string, bool>> spectrahedra;

    Point objFunction;

    SdpaFormatManager<NT> sdpaFormatManager;

    std::ifstream in1;
    SPECTRAHEDRON spectrahedron1;
    in1.open("spectra_data/sdp_prob_200_15.txt", std::ifstream::in);
    sdpaFormatManager.loadSDPAFormatFile(in1, spectrahedron1, objFunction);
    spectrahedra.push_back(std::make_tuple(spectrahedron1, "S_200_15", true));

    std::ifstream in2;
    SPECTRAHEDRON spectrahedron2;
    in2.open("spectra_data/sdp_prob_400_20.txt", std::ifstream::in);
    sdpaFormatManager.loadSDPAFormatFile(in2, spectrahedron2, objFunction);
    spectrahedra.push_back(std::make_tuple(spectrahedron2, "S_400_20", true));

    std::ifstream in3;
    SPECTRAHEDRON spectrahedron3;
    in3.open("spectra_data/sdp_prob_600_25.txt", std::ifstream::in);
    sdpaFormatManager.loadSDPAFormatFile(in3, spectrahedron3, objFunction);
    spectrahedra.push_back(std::make_tuple(spectrahedron3, "S_600_25", true));

    std::ifstream in4;
    SPECTRAHEDRON spectrahedron4;
    in4.open("spectra_data/sdp_prob_800_30.txt", std::ifstream::in);
    sdpaFormatManager.loadSDPAFormatFile(in4, spectrahedron4, objFunction);
    spectrahedra.push_back(std::make_tuple(spectrahedron4, "S_800_30", true));

    std::ifstream in5;
    SPECTRAHEDRON spectrahedron5;
    in5.open("spectra_data/sdp_prob_1000_35.txt", std::ifstream::in);
    sdpaFormatManager.loadSDPAFormatFile(in5, spectrahedron5, objFunction);
    spectrahedra.push_back(std::make_tuple(spectrahedron5, "S_1000_35", true));


    SPECTRAHEDRON P;
    std::string name;
    std::ofstream outfile;
    NT step_size = 0;
    std::pair<Point, NT> inner_ball;

    for (std::tuple<SPECTRAHEDRON, std::string, bool> spectrahedron_tuple : spectrahedra) {
        P = std::get<0>(spectrahedron_tuple);
        name = std::get<1>(spectrahedron_tuple);
        std::cout << name << std::endl;
        outfile.open("results_" + name + "_new.txt");
        //P.normalize();
        //inner_ball = P.ComputeInnerBall();
        step_size = 0.5;
        for (unsigned int walk_length = P.dimension()/2; walk_length <= P.dimension(); walk_length += P.dimension() / 2) {
            results = benchmark_spectrahedron_sampling<NT, SPECTRAHEDRON>(P, step_size, walk_length, std::get<2>(spectrahedron_tuple));
            outfile << results[0];
            outfile << results[1];
        }
        outfile.close();
    }

}

template <typename NT>
void call_test_optimization() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;

    Hpolytope P = generate_cube<Hpolytope>(100, false);

    Point coeffs = Point::all_ones(100);

    benchmark_polytope_linear_program_optimization<NT, Hpolytope>(coeffs, P);

}

template <typename NT>
void call_test_benchmark_convex_body() {
  typedef Cartesian<NT>    Kernel;
  typedef typename Kernel::Point    Point;
  typedef ConvexBody<Point> Convexbody;
  typedef boost::mt19937 RNGType;
  typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
  typedef InnerBallFunctor::GradientFunctor<Point> NegativeGradientFunctor;
  typedef InnerBallFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
  typedef GeneralizedLeapfrogODESolver<Point, NT, Convexbody, NegativeGradientFunctor> Solver;
  typedef typename Convexbody::MT MT;
  typedef typename Convexbody::VT VT;

  unsigned int dim = 1000;
  // Convexbody P = generate_unit_ball<Convexbody>(dim);
  Convexbody P = generate_unit_ball_intersect_logsumexp<Convexbody>(dim);

  unsigned int max_draws = 40000;
  unsigned int num_burns = max_draws / 3;

  Point x0(dim);
  NT R0 = NT(dim);

  unsigned int walk_length = 1000;

  SimulationStats<NT> hmc_stats;

  // Random number generator
  RandomNumberGenerator rng(1);

  // Declare oracles
  InnerBallFunctor::parameters<NT, Point> params(x0, R0);

  NegativeGradientFunctor F(params);
  NegativeLogprobFunctor f(params);

  int max_actual_draws = max_draws - num_burns;
  unsigned int min_ess = 0;

  MT samples;
  samples.resize(dim, max_actual_draws);
  NT ETA;
  NT max_psrf;

  std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;

  std::cerr << std::endl;
  std::cerr << "HMC Sampling" << std::endl;


  HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);

  HamiltonianMonteCarloWalk::Walk
    <Point, Convexbody, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
    hmc(&P, x0, F, f, hmc_params);

  min_ess = 0;

  std::cerr << "Hamiltonian Monte Carlo (Gaussian Density)" << std::endl;

  hmc.solver->eta = 0.001;

  std::cerr << "Burn-in" << std::endl;

  for (unsigned int i = 0; i < num_burns; i++) {
    if (i % 1000 == 0) std::cerr << ".";
    hmc.apply(rng, walk_length);
  }

  hmc.disable_adaptive();
  std::cerr << std::endl;
  std::cerr << "Sampling" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  for (unsigned int i = 0; i < max_actual_draws; i++) {
    hmc.apply(rng, walk_length);
    std::cout << hmc.x.getCoefficients().transpose() << std::endl;

    samples.col(i) = hmc.x.getCoefficients();
    if (i % 1000 == 0 && i > 0) std::cerr << ".";
  }
  stop = std::chrono::high_resolution_clock::now();

  ETA = (NT) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

  std::cerr << std::endl;
  print_diagnostics<NT, VT, MT>(samples, min_ess, std::cerr);
  std::cerr << "Average time per sample: " << ETA / max_actual_draws << "us" << std::endl;
  std::cerr << "Average time per independent sample: " << ETA / min_ess << "us" << std::endl;
  std::cerr << "Average number of reflections: " <<
      (1.0 * hmc.solver->num_reflections) / hmc.solver->num_steps << std::endl;
  std::cerr << "Step size (final): " << hmc.solver->eta << std::endl;
  std::cerr << "Discard Ratio: " << hmc.discard_ratio << std::endl;
  std::cerr << "Average Acceptance Probability: " << exp(hmc.average_acceptance_log_prob) << std::endl;
  std::cerr << std::endl;
  std::cerr << "Min ESS" << min_ess << std::endl;

  max_psrf = check_interval_psrf<NT, VT, MT>(samples);

  hmc_stats.method = "HMC";
  hmc_stats.walk_length = walk_length;
  hmc_stats.min_ess = min_ess;
  hmc_stats.max_psrf = max_psrf;
  hmc_stats.time_per_draw = ETA / max_actual_draws;
  hmc_stats.time_per_independent_sample = ETA / min_ess;
  hmc_stats.average_number_of_reflections = (1.0 * hmc.solver->num_reflections) / hmc.solver->num_steps;
  hmc_stats.step_size = hmc.solver->eta;
  hmc_stats.average_acceptance_log_prob  = exp(hmc.average_acceptance_log_prob);

}

TEST_CASE("hmc") {
    call_test_hmc<double>();
}

TEST_CASE("uld") {
    call_test_uld<double>();
}

TEST_CASE("optimization") {
    call_test_optimization<double>();
}

TEST_CASE("benchmark_hmc") {
    call_test_benchmark_hmc<double>(false);
}

TEST_CASE("benchmark_hmc_truncated") {
    call_test_benchmark_hmc<double>(true);
}

TEST_CASE("benchmark_polytopes_grid_search") {
    call_test_benchmark_polytopes_grid_search<double>();
}

TEST_CASE("benchmark_convex_body") {
    call_test_benchmark_convex_body<double>();
}

TEST_CASE("benchmark_spectrahedra_grid_search") {
    call_test_benchmark_spectrahedra_grid_search<double>();
}
