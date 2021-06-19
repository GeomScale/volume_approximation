// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

//Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2018 and 2019 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "sampling/sampling.hpp"
#include "ode_solvers/ode_solvers.hpp"
#include "oracle_functors_rcpp.h"

template <
  typename Solver,
  typename MT
>
void run_ode_solver(
  Solver &solver,
  unsigned int &order,
  unsigned int &num_steps,
  unsigned int &dimension,
  Rcpp::List &results) {

  std::vector<MT> results_temp;

  for (unsigned int i = 0; i < order; i++) {
    MT temp_result;
    temp_result.resize(dimension, num_steps);
    results_temp.push_back(temp_result);
  }

  for (unsigned int i = 0; i < num_steps; i++) {
    for (unsigned int j = 0; j < order; j++) {
      results_temp[j].col(i) = solver.xs[j].getCoefficients();
    }
    solver.step(i, true);
  }

  for (unsigned int i = 0; i < order; i++) {
    std::ostringstream stringStream;
    stringStream << "x_" << i + 1;
    std::string state_name = stringStream.str();

    results.push_back(Rcpp::wrap(results_temp[i]), state_name.c_str());

  }

}

//' Solve an ODE of the form dx^n / dt^n = F(x, t)
//'
//' @param n The number of steps.
//' @param step_size The step size.
//' @param order The ODE order (default is n = 1)
//' @param dimension The dimension of each derivative
//' @param initial_time The initial time
//' @param F The function oracle F(x, t) in the ODE.
//' @param method The method to be used
//' @param initial_conditions The initial conditions provided to the solver. Must be provided in a list with keys "x_1", ..., "x_n" and column vectors as values. The state "x_n" represents the (n-1)-th order derivative with respect to time
//' @param domains A list of n H-polytopes with keys "P_1", "P_2", ..., "P_n" that correspond to each derivative's domain
//'
//' @return A list which contains elements "x_1", ..., "x_n" representing each derivative results. Each "x_i" corresponds to a d x n matrix where each column represents a certain timestep of the solver.
//'
//' @examples
//' # Please visit the examples directory on examples demonstrating usage of the ODE solvers.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List ode_solve(Rcpp::Nullable<unsigned int> n,
                     Rcpp::Nullable<double> step_size,
                     Rcpp::Nullable<unsigned int> order,
                     Rcpp::Nullable<unsigned int> dimension,
                     Rcpp::Nullable<double> initial_time,
                     Rcpp::Function F,
                     Rcpp::Nullable<Rcpp::CharacterVector> method,
                     Rcpp::Nullable<Rcpp::List> domains = R_NilValue,
                     Rcpp::Nullable<Rcpp::List> initial_conditions = R_NilValue)
{

  typedef double NT;
  typedef Cartesian<NT>    Kernel;
  typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
  typedef typename Kernel::Point    Point;
  typedef HPolytope <Point> Hpolytope;
  typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
  typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
  typedef std::vector<Point> pts;
  typedef std::vector<Hpolytope*> hp_bounds;
  typedef RcppFunctor::GradientFunctor<Point> func;

  unsigned int n_ = Rcpp::as<unsigned int>(n);
  NT step_size_ = Rcpp::as<NT>(step_size);
  unsigned int order_ = Rcpp::as<unsigned int>(order);
  std::string method_ = Rcpp::as<std::string>(method);
  unsigned int dim = Rcpp::as<unsigned int>(dimension);
  NT initial_time_ = Rcpp::as<NT>(initial_time);

  // Create functors
  RcppFunctor::parameters<NT> rcpp_functor_params(1, 1, order_);

  // Create C++ functor
  func F_(rcpp_functor_params, F, false);

  // Initialize initial conditions
  pts initial_conditions_;
  VT temp_initial_condition;

  for (unsigned int i = 1; i <= order_; i++) {
    std::ostringstream stringStream;
    stringStream << "x_" << i;
    std::string state_name = stringStream.str();

    if (Rcpp::as<Rcpp::List>(initial_conditions).containsElementNamed(state_name.c_str())) {
      temp_initial_condition = Rcpp::as<VT>(Rcpp::as<Rcpp::List>(initial_conditions)[state_name.c_str()]);
      Point p(temp_initial_condition);
      initial_conditions_.push_back(p);
    } else {
      Point p(dim);
      initial_conditions_.push_back(p);
    }

  }

  unsigned int i = 0;
  unsigned int t = 0;
  F_(i, initial_conditions_, t);

  hp_bounds domains_;

  // Initialize domains
  // TODO Test and add other polytope types
  for (unsigned int i = 0; i < order_; i++) {
    std::ostringstream stringStream;
    stringStream << "P_" << i + 1;
    std::string domain_name = stringStream.str();

    if (Rcpp::as<Rcpp::List>(domains).containsElementNamed(domain_name.c_str())) {

      Hpolytope HP(dim, Rcpp::as<MT>(
                    Rcpp::as<Rcpp::Reference>(Rcpp::as<Rcpp::List>(domains)
                    [domain_name.c_str()]).field("A")),
                   Rcpp::as<VT>(
                     Rcpp::as<Rcpp::Reference>(Rcpp::as<Rcpp::List>(domains)
                     [domain_name.c_str()]).field("b"))
              );

      HP.normalize();

      Hpolytope *HP_ref = &HP;

      if (!HP_ref->is_in(initial_conditions_[i])) {
        std::ostringstream errStream;
        errStream << "Initial condition out of bounds for state " << i + 1 << ".\n";
        throw Rcpp::exception(errStream.str().c_str());
      }
      domains_.push_back(HP_ref);
    } else {
      domains_.push_back(NULL);
    }
  }

  Rcpp::List results;

  if (method_ == "euler") {
    EulerODESolver<Point, NT, Hpolytope, func> euler_solver =
      EulerODESolver<Point, NT, Hpolytope, func>
      (initial_time_, step_size_, initial_conditions_, F_, domains_);
    run_ode_solver<EulerODESolver<Point, NT, Hpolytope, func>, MT>(euler_solver, order_, n_, dim, results);
  } else if (method_ == "leapfrog") {
    if (order_ % 2 == 1) {
      throw Rcpp::exception("Leapfrog is an even order solver.");
    }
    LeapfrogODESolver<Point, NT, Hpolytope, func> leapfrog_solver =
      LeapfrogODESolver<Point, NT, Hpolytope, func>
      (initial_time_, step_size_, initial_conditions_, F_, domains_);
    run_ode_solver<LeapfrogODESolver<Point, NT, Hpolytope, func>, MT>(leapfrog_solver, order_, n_, dim, results);
  } else if (method_ == "runge_kutta") {
    RKODESolver<Point, NT, Hpolytope, func> rk_solver =
      RKODESolver<Point, NT, Hpolytope, func>
      (initial_time_, step_size_, initial_conditions_, F_, domains_);
    run_ode_solver<RKODESolver<Point, NT, Hpolytope, func>, MT>(rk_solver, order_, n_, dim, results);
  } else if (method_ == "richardson") {
    RichardsonExtrapolationODESolver<Point, NT, Hpolytope, func> r_solver =
      RichardsonExtrapolationODESolver<Point, NT, Hpolytope, func>
      (initial_time_, step_size_, initial_conditions_, F_, domains_);
    run_ode_solver<RichardsonExtrapolationODESolver<Point, NT, Hpolytope, func>, MT>(r_solver, order_, n_, dim, results);
  } else {
    throw Rcpp::exception("Unrecognized solver. Aborting.");
  }

  return results;
}
