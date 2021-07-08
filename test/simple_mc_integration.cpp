// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Suraj Choubey, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include "simple_MC_integration.hpp"
#include "Eigen/Eigen"
#include <vector>
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/known_polytope_generators.h"
#include "ode_solvers/oracle_functors.hpp"
#include "random_walks/random_walks.hpp"
#include <iostream>
#include <fstream>
#include "misc.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef std::vector<Point> Points;
typedef HPolytope<Point> HPOLYTOPE;
typedef VPolytope<Point> VPOLYTOPE;
typedef boost::mt19937 RNGType;
typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;

NT exp_normsq(Point X) {
	return exp(-X.squared_length()) ;
}

NT simple_polynomial_1D(Point X) {
	return (X[0] - 1) * (X[0] - 2) * (X[0] - 3);
}

NT logx_natural_1D(Point X) {
	return log(X[0]);
}

NT rooted_squaresum(Point X) {
	return sqrt(X.squared_length());
}

NT one_sqsum(Point X) {
	return 1 - X.squared_length();
}

template <typename NT>
void test_values (NT computed, NT expected, NT exact) {
	std::cout << "Computed integration value = " << computed << std::endl;
	std::cout << "Expected integration value = " << expected << std::endl;
	std::cout << "Exact integration value = " << exact << std::endl;
	std::cout << "Relative error (expected) = " << std::abs((computed - expected)/expected) << std::endl;
	std::cout << "Relative error (exact) = " << std::abs((computed - exact)/exact) << std::endl ;
    CHECK(((std::abs((computed - expected)/expected) < 0.00001) || (std::abs((computed - exact)/exact) < 0.2)));
}

template <typename NT>
void call_test_simple_mc_integration_over_rectangles() {
	
	NT integration_value;
	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER RECTANGLES USING UNIFORM RANDOM WALKS\n";

	Limit LL{-1};
	Limit UL{6};
	integration_value = simple_mc_integrate <AcceleratedBilliardWalk> (simple_polynomial_1D, 1, 100000, CB, LL, UL);
	test_values(integration_value, 39.7, 40.25);

	Limit LL1{0.5};
	Limit UL1{10};
	integration_value = simple_mc_integrate <BilliardWalk> (logx_natural_1D, 1, 100000, CB, LL1, UL1);
	test_values(integration_value, 13.65, 13.872);

	Limit LL2{-1, -1};
	Limit UL2{1, 1};
	integration_value = simple_mc_integrate <AcceleratedBilliardWalk> (rooted_squaresum, 2, 100000, SOB, LL2, UL2);
	test_values(integration_value, 2.99, 3.0607);

	integration_value = simple_mc_integrate <BilliardWalk> (exp_normsq, 5, 100000, SOB);
	test_values(integration_value, 7.49, 7.46);

	integration_value = simple_mc_integrate <AcceleratedBilliardWalk> (exp_normsq, 8, 100000, SOB);
	test_values(integration_value, 24.8, 24.76);

	integration_value = simple_mc_integrate <BilliardWalk> (exp_normsq, 10, 100000, SOB);
	test_values(integration_value, 54.8, 55.25);

}

template <typename NT>
void call_test_simple_mc_integration_over_cubes() {

	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER CUBES USING UNIFORM RANDOM WALKS\n";

	NT integration_value;
	HPOLYTOPE HP;

	HP = generate_cube <HPOLYTOPE> (2, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (exp_normsq, HP, 100000, SOB, 10, 0.01);
	test_values(integration_value, 2.20, 2.230);

	// For 2D Polytope shifted to (1,1) from origin
	std::vector<NT> Origin{1, 1};
	Point newOrigin(2, Origin);
	integration_value = simple_mc_polytope_integrate <AcceleratedBilliardWalk, HPOLYTOPE> (exp_normsq, HP, 100000, SOB, 1, 0.01, newOrigin);
	test_values(integration_value, 0.78, 0.777);

	HP = generate_cube <HPOLYTOPE> (10, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (exp_normsq, HP, 100000, SOB);
	test_values(integration_value, 54.7, 55.25);

	HP = generate_cube <HPOLYTOPE> (15, false);
	integration_value = simple_mc_polytope_integrate <AcceleratedBilliardWalk, HPOLYTOPE> (exp_normsq, HP, 100000, SOB);
	test_values(integration_value, 405.9, 410.690);

	HP = generate_cube <HPOLYTOPE> (20, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (exp_normsq, HP, 100000, SOB);
	test_values(integration_value, 3050.0, 3052.71);

	// Reading a H-Polytope from *.ine file for 20 Dimensions
	// std::string fileName("cube10.ine");
	// std::ifstream inp;
	// std::vector<std::vector<NT>> Pin;
	// inp.open(fileName, std::ifstream::in);
	// read_pointset(inp,Pin);
	// HPOLYTOPE HP4(Pin);
	// integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (exp_normsq, 2, HP4, 100000, SOB);
	// std::cout << "Integration value: " << integration_value << std::endl;
	// test_values(integration_value, expected, exact);
	// inp.close();
}

template <typename NT>
void call_test_simple_mc_integration_over_simplices() {

	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER SIMPLICES USING UNIFORM RANDOM WALKS\n";

	NT integration_value;
	HPOLYTOPE HP;

	HP = generate_simplex <HPOLYTOPE> (1, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, CB, 10, 0.01);
	test_values(integration_value, 0.67, 0.666);
	
	HP = generate_simplex <HPOLYTOPE> (2, false);
	integration_value = simple_mc_polytope_integrate <AcceleratedBilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB, 10, 0.01);
	test_values(integration_value, 0.34, 0.333);

	HP = generate_simplex <HPOLYTOPE> (3, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB);
	test_values(integration_value, 0.116, 0.1166);

	HP = generate_simplex <HPOLYTOPE> (5, false);
	integration_value = simple_mc_polytope_integrate <AcceleratedBilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB);
	test_values(integration_value, 0.00656, 0.0063492);

	HP = generate_simplex <HPOLYTOPE> (7, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB);
	test_values(integration_value, 0.000159, 0.000159832);

}

template <typename NT>
void call_test_simple_mc_integration_over_product_simplices() {

	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER PRODUCT SIMPLICES USING UNIFORM RANDOM WALKS\n";

	NT integration_value;
	HPOLYTOPE HP;

	HP = generate_prod_simplex <HPOLYTOPE> (1, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, CB, 10, 0.01);
	test_values(integration_value, 0.334, 0.333);
	
	HP = generate_prod_simplex <HPOLYTOPE> (2, false);
	integration_value = simple_mc_polytope_integrate <AcceleratedBilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB);
	test_values(integration_value, 0.0834, 0.0833);

	HP = generate_prod_simplex <HPOLYTOPE> (3, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB);
	test_values(integration_value, 0.0110, 0.01111);

	HP = generate_prod_simplex <HPOLYTOPE> (5, false);
	integration_value = simple_mc_polytope_integrate <AcceleratedBilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB);
	test_values(integration_value, 0.36e-4, 0.36375e-4);

	HP = generate_prod_simplex <HPOLYTOPE> (7, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (one_sqsum,  HP, 100000, SOB);
	test_values(integration_value, 0.235e-7, 0.24079e-7);

}

template <typename NT>
void call_test_simple_mc_integration_over_cross_polytopes() {

	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER CROSS POLYTOPES USING UNIFORM RANDOM WALKS\n";

	NT integration_value;
	HPOLYTOPE HP;

	HP = generate_cross <HPOLYTOPE> (1, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, CB, 10, 0.01);
	test_values(integration_value, 1.334, 1.333333);
	
	HP = generate_cross <HPOLYTOPE> (2, false);
	integration_value = simple_mc_polytope_integrate <AcceleratedBilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB, 10, 0.01);
	test_values(integration_value, 1.334, 1.33333);

	HP = generate_cross <HPOLYTOPE> (3, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB);
	test_values(integration_value, 0.935000, 0.933333);

	HP = generate_cross <HPOLYTOPE> (5, false);
	integration_value = simple_mc_polytope_integrate <AcceleratedBilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB);
	test_values(integration_value, 0.200000, 0.203174);

	HP = generate_cross <HPOLYTOPE> (7, false);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB);
	test_values(integration_value, 0.020000, 0.020458);

}

template <typename NT>
void call_test_simple_mc_integration_over_birkhoff_polytopes() {

	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER BIRKHOFF POLYTOPES USING UNIFORM RANDOM WALKS\n";

	NT integration_value;
	HPOLYTOPE HP;

	HP = generate_birkhoff <HPOLYTOPE> (2);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, CB, 10, 0.01);
	test_values(integration_value, 0.67, 0.6666);

	HP = generate_birkhoff <HPOLYTOPE> (3);
	integration_value = simple_mc_polytope_integrate <AcceleratedBilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB);
	test_values(integration_value, 0.0470, 0.04722);
	
	HP = generate_birkhoff <HPOLYTOPE> (4);
	integration_value = simple_mc_polytope_integrate <BilliardWalk, HPOLYTOPE> (one_sqsum, HP, 100000, SOB);
	test_values(integration_value, 0.000150, 0.000164);

}

TEST_CASE("rectangle") {
    call_test_simple_mc_integration_over_rectangles<double>();
}

TEST_CASE("cube") {
    call_test_simple_mc_integration_over_cubes<double>();
}

TEST_CASE("simplex") {
    call_test_simple_mc_integration_over_simplices<double>();
}

TEST_CASE("prod_simplex") {
	call_test_simple_mc_integration_over_product_simplices<double>();
}

TEST_CASE("cross") {
	call_test_simple_mc_integration_over_cross_polytopes<double>();
}

TEST_CASE("birkhoff") {
	call_test_simple_mc_integration_over_birkhoff_polytopes<double>();
}
