// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "preprocess/min_sampling_covering_ellipsoid_rounding.hpp"
#include "preprocess/svd_rounding.hpp"
#include "preprocess/max_inscribed_ellipsoid_rounding.hpp"
#include "extractMatPoly.h"

template 
<
    typename MT,
    typename VT, 
    typename WalkType, 
    typename Polytope, 
    typename Point, 
    typename NT, 
    typename RNGType
>
std::tuple<MT, VT, NT> apply_rounding(Polytope &P,
                                      std::string const& method_rcpp,
                                      unsigned int const& walkL,
                                      std::pair<Point, NT> &InnerBall, 
                                      RNGType &rng) 
{
    std::tuple<MT, VT, NT> round_res;
    if (method_rcpp.compare(std::string("min_ellipsoid")) == 0) {
        round_res = min_sampling_covering_ellipsoid_rounding<WalkType, MT, VT>(P, InnerBall, walkL, rng);
    } else if (method_rcpp.compare(std::string("svd")) == 0) {
        round_res = svd_rounding<WalkType, MT, VT>(P, InnerBall, walkL, rng);
    } else {
        throw Rcpp::exception("Unknown method!");
    }
    return round_res;
}

//' Internal rcpp function for the rounding of a convex polytope
//'
//' @param P A convex polytope (H- or V-representation or zonotope).
//' @param method Optional. The method to use for rounding, a) \code{'min_ellipsoid'} for the method based on mimimmum volume enclosing ellipsoid of a uniform sample from P, b) \code{'max_ellipsoid'} for the method based on maximum volume enclosed ellipsoid in P, (c) \code{'svd'} for the method based on svd decomposition. The default method is \code{'min_ellipsoid'} for all the representations.
//' @param seed Optional. A fixed seed for the number generator.
//'
//' @keywords internal
//'
//' @return A numerical matrix that describes the rounded polytope, a numerical matrix of the inverse linear transofmation that is applied on the input polytope, the numerical vector the the input polytope is shifted and the determinant of the matrix of the linear transformation that is applied on the input polytope.
// [[Rcpp::export]]
Rcpp::List rounding (Rcpp::Reference P, 
                     Rcpp::Nullable<std::string> method = R_NilValue,
                     Rcpp::Nullable<double> seed = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;
    typedef Zonotope<Point> zonotope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    unsigned int n = P.field("dimension"), walkL = 2, type = P.field("type");
    std::string method_rcpp = std::string("min_ellipsoid");
    if(method.isNotNull()) {
        method_rcpp =  Rcpp::as<std::string>(method);
        if (method_rcpp.compare(std::string("max_ellipsoid")) == 0 && type != 1) {
            Rcpp::exception("This method can not be used for V- or Z-polytopes!");
        }
    }

    RNGType rng(n);
    if (seed.isNotNull()) {
        unsigned seed_rcpp = Rcpp::as<double>(seed);
        rng.set_seed(seed_rcpp);
    }

    std::pair <Point, NT> InnerBall;
    Rcpp::NumericMatrix Mat;

    std::tuple<MT, VT, NT> round_res;
    switch (type) {
        case 1: {
            // Hpolytope

            if (Rcpp::as<MT>(P.field("Aeq")).rows() > 0) {
                throw Rcpp::exception("volesti supports rounding for full dimensional polytopes. Maybe call function get_full_dimensional_polytope()");
            } 
            Hpolytope HP(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            InnerBall = HP.ComputeInnerBall();
            if (InnerBall.second < 0.0) throw Rcpp::exception("Unable to compute a feasible point.");
            if (method_rcpp.compare(std::string("max_ellipsoid")) == 0) {
                round_res = max_inscribed_ellipsoid_rounding<MT, VT, NT>(HP, InnerBall.first);
            } else {
                round_res = apply_rounding<MT, VT, AcceleratedBilliardWalk>(HP, method_rcpp, walkL, InnerBall, rng);
            }
            Mat = extractMatPoly(HP);
            break;
        }
        case 2: {
            // Vpolytope
            Vpolytope VP(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            InnerBall = VP.ComputeInnerBall();
            if (InnerBall.second < 0.0) throw Rcpp::exception("Unable to compute a feasible point.");
            round_res = apply_rounding<MT, VT, BilliardWalk>(VP, method_rcpp, walkL, InnerBall, rng);
            Mat = extractMatPoly(VP);
            break;
        }
        case 3: {
            // Zonotope
            zonotope ZP(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            InnerBall = ZP.ComputeInnerBall();
            if (InnerBall.second < 0.0) throw Rcpp::exception("Unable to compute a feasible point.");
            round_res = apply_rounding<MT, VT, BilliardWalk>(ZP, method_rcpp, walkL, InnerBall, rng);
            Mat = extractMatPoly(ZP);
            break;
        }
        case 4: {
            throw Rcpp::exception("volesti does not support rounding for this representation currently.");
        }
    }

    return Rcpp::List::create(Rcpp::Named("Mat") = Mat, Rcpp::Named("T") = Rcpp::wrap(std::get<0>(round_res)),
                              Rcpp::Named("shift") = Rcpp::wrap(std::get<1>(round_res)),
                              Rcpp::Named("round_value") = std::get<2>(round_res));

}
