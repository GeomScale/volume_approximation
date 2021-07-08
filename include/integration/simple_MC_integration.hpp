// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Suraj Choubey, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

// Monte Carlo Integration algorithm used here : https://en.wikipedia.org/wiki/Monte_Carlo_integration#Overview

#ifndef SIMPLE_MC_INTEGRATION_HPP
#define SIMPLE_MC_INTEGRATION_HPP

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "convex_bodies/hpolytope.h"
#include "Eigen/Eigen"
#include "generators/known_polytope_generators.h"
#include "boost_random_number_generator.hpp"
#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "misc.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef std::vector<Point> Points;
typedef HPolytope<Point> HPOLYTOPE;
typedef boost::mt19937 RNGType;
typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
typedef typename HPolytope<Point>::MT MT;
typedef typename HPolytope<Point>::VT VT; 

typedef const unsigned int Uint;  // Positive constant value for no of samples & dimensions
enum volumetype {CB ,CG ,SOB}; // Volume type for polytope
typedef typename std::vector<NT> Limit; // Standard way for user to use limits 
// E.g. Limits LL{0.5, 1.5, 2.5} ; Limits UL{1.2, 1.8 , 2.8 } ;

const Limit lt{0}; // To initialize non-initialized limits
const Point pt{0}; // To initialize non-initialized points

// To check if two n-dimensional points ensure valid limits in integration
template
<
    typename Point = Point,
    typename NT = NT
>
bool valid_limits(Point LL, Point UL) {
    if (UL.dimension() == LL.dimension()) {
        for (int i = 0; i<LL.dimension(); i++) {
            if (LL[i] > UL[i]) {
                std::cerr << "Invalid integration limits\n";
                return false;
            }
        }
        return true;
    } else {
        std::cerr << "Invalid integration limits\n";
        return false;
    }
}

// Initialize Limit Point
template
<
    typename Point = Point,
    typename NT = NT
>
Point init_limit(Limit L, Uint dim) {
    Point pt(dim);
    for (int i=0; i<dim; i++) {
        pt.set_coord(i, L[i]);
    }
    return pt;
}

//Initialize to [-1,1]^n
template
<
    typename Point = Point,
    typename NT = NT
>
void initiate_unit_limits(Point& LL, Point& UL, int dim) {
    LL.set_dimension(dim);
    UL.set_dimension(dim);
    LL.set_to_origin();
    UL.set_to_origin();

    for (int i=0 ; i<dim; i++) {
        LL.set_coord(i, -1);
        UL.set_coord(i, 1);
    }
}

// Simple MC Integration Over Polytopes
template 
<
    typename WalkType = BallWalk,
    typename Polytope = HPOLYTOPE,
    typename RNG = RandomNumberGenerator,
    typename NT = NT,
    typename Functor
>
NT simple_mc_polytope_integrate(Functor Fx, 
                                Polytope &P, 
                                Uint N = 10000, 
                                volumetype voltype = SOB, 
                                int walk_length = 1, 
                                NT e = 0.1, 
                                Point Origin = pt) 
{

    Uint dim = P.dimension();
    // P.print();

    // Check if ShiftPoint is shifted with accurate dimensions
    if (Origin.dimension() == 0 && dim > 0) {
        Origin.set_dimension(dim); 
        Origin.set_to_origin();
    } else if (Origin.dimension() != dim && dim > 0) {
        std::cerr << "Polytope Dimension != Shiftpoint Dimension" << std::endl;
        return -1;
    }

    // std::cout << "Origin.dimension() = " << Origin.dimension() << std::endl;
    // std::cout << "P.dimension() = " << P.dimension() << std::endl;

    // Volume calculation for HPolytope
    NT volume;
    
    switch (voltype) {
    case CB:     
        volume = volume_cooling_balls <BallWalk, RNG, Polytope> (P, e, walk_length).second; 
        break;
    case CG: 
        volume = volume_cooling_gaussians <GaussianBallWalk, RNG, Polytope> (P, e, walk_length);
        break;
    case SOB: 
        volume = volume_sequence_of_balls <BallWalk, RNG, Polytope> (P, e, walk_length);
        break;
    default:
        std::cerr << "Error in volume type: CB / SOB / CG" << std::endl;
        return -1;
    }

    // std::cout << "Volume of the convex body = " << volume << std::endl;

    // For implementing Uniform Walks
    RNG rng(1);
    std::pair <Point, NT> inner_ball = P.ComputeInnerBall();
    Point x0 = inner_ball.first;
    typename WalkType::template Walk<Polytope, RNG> walk(P, x0, rng);

    // For storing sampled points
    Points points; 
    NT sum = 0;

    // Applying and walking through Uniform Walks + Storing Points in Vector<Point>
    for (int i = 0; i < N; i++) {
        walk.apply(P, x0, walk_length, rng);
        sum += Fx(x0 + Origin);
        
        // (x0 + Origin).print();
    } 

    // Integration Value
    NT integration_value = volume * sum / N ;
    return integration_value;
}

// Simple MC Integration over Hyper-Rectangles
template
<
    typename WalkType = BallWalk,
    typename RNG = RandomNumberGenerator,
    typename NT = NT,
    typename Functor
>
NT simple_mc_integrate (Functor Fx, 
                        Uint dim, 
                        Uint N = 10000, 
                        volumetype voltype = SOB, 
                        Limit LowLimit = lt, 
                        Limit UpLimit = lt, 
                        int walk_length = 10, 
                        NT e = 0.1) 
{

    // Setting up integration limits
    Point LL, UL;
    if (LowLimit.size() == 1 && UpLimit.size() == 1 && LowLimit[0] == 0 && UpLimit[0] == 0) { 
        initiate_unit_limits(LL, UL, dim); 
    } else if (LowLimit.size() == UpLimit.size() && LowLimit.size() == dim) { 
        LL = init_limit <Point, NT> (LowLimit, dim); 
        UL = init_limit <Point, NT> (UpLimit, dim); 
    } else { 
        std::cerr << "Invalid limits entered"; 
        return -1; 
    }
 
    NT sum = 0;

    if (valid_limits(LL, UL)) {

        // Creating an MT & VT for HPolytope(Hyper-Rectangle) for integration limits using LL & UL
        MT mt(dim*2, dim); 
        mt.setZero();
        VT vt(dim*2); 
        vt.setZero();

        for (int i=0; i<dim; i++) {
            mt(i,i) = 1;
            vt(i) = UL[i];
            mt(dim+i,i) = -1;
            vt(dim+i) = -1 * LL[i];
        }

        // Initialization of H-Polytope
        HPOLYTOPE P(dim, mt, vt);
        // P.print();

        NT integration_value = simple_mc_polytope_integrate <WalkType, HPOLYTOPE> (Fx, P, N, voltype, walk_length, e);
        return integration_value;

    } else {
        std::cerr << "Invalid integration limits" << std::endl;
        return -1;
    }
}

#endif

// Hyper-Rectangle volume calculator in n-dimensions
// template
// <
//     typename Point = Point,
//     typename NT = NT
// >
// NT hyper_rect_volume (Point LL, Point UL) {

//     NT product = 1;
//     if( valid_limits(LL, UL) ) {
//         for(int i=0; i<LL.dimension() ; ++i) {
//             product = product * (UL[i] - LL[i]);
//         }
//         return product;
//     }
//     else return -1;
// }

// To sample a point between two n-dimensional points using inbuilt random sampling
// Please use `srand(time(0));` while invoking this function
// Point sample_point_bw_limits (Point LL, Point UL) {
    
//     Point sample_point(LL.dimension());
//     for (int i=0; i<LL.dimension() ; ++i) {
//         sample_point.set_coord(i, LL[i] + (NT)(rand()) / ((NT)(RAND_MAX/(UL[i] - LL[i]))));
//     }
//     return sample_point;
// }

// Setting up params for random walks
// std::pair <Point, NT> inner_ball = P.ComputeInnerBall();
// RNG rng(1);
// Point x0 = inner_ball.first;
// typename WalkType::template Walk <HPOLYTOPE, RNG> walk(P, x0, rng);

// for (int i = 0; i<=N; i++) {
//     walk.apply(P, x0, walk_length, rng);
//     sum += Fx(x0);
// }    

// NT volume = hyper_rect_volume(LL,UL);
// std::cout << volume << std::endl;

// NT integration_value = volume * sum / N ;
// return integration_value;
