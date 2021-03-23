// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/rotating.hpp"
#include "extractMatPoly.h"

//'  An internal Rccp function for the random rotation of a convex polytope
//'
//' @param P A convex polytope (H-, V-polytope or a zonotope).
//' @param T Optional. A rotation matrix.
//' @param seed Optional. A fixed seed for the random linear map generator.
//'
//' @keywords internal
//'
//' @return A matrix that describes the rotated polytope
// [[Rcpp::export]]
Rcpp::NumericMatrix rotating (Rcpp::Reference P, Rcpp::Nullable<Rcpp::NumericMatrix> T = R_NilValue,
                              Rcpp::Nullable<int> seed = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;
    typedef Zonotope<Point> zonotope;
    typedef IntersectionOfVpoly<Vpolytope, RNGType> InterVP;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    MT TransorfMat;
    Rcpp::NumericMatrix Mat;
    unsigned int n = P.field("dimension"), type = P.field("type");

    int seed_rcpp = (!seed.isNotNull()) ? std::chrono::system_clock::now()
                                      .time_since_epoch().count()
                                    : Rcpp::as<int>(seed);

    switch (type) {
        case 1: {
            // Hpolytope
            Hpolytope HP(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            if (T.isNotNull()) {
                TransorfMat = Rcpp::as<MT>(T);
                HP.linear_transformIt(TransorfMat.inverse());
            } else {
                TransorfMat = rotating < MT > (HP, seed_rcpp);
            }
            Mat = extractMatPoly(HP);
            break;
        }
        case 2: {
            // Vpolytope
            Vpolytope VP(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            if (T.isNotNull()) {
                TransorfMat = Rcpp::as<MT>(T);
                VP.linear_transformIt(TransorfMat.inverse());
            } else {
                TransorfMat = rotating < MT > (VP, seed_rcpp);
            }
            Mat = extractMatPoly(VP);
            break;
        }
        case 3: {
            // Zonotope
            zonotope ZP(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            if (T.isNotNull()) {
                TransorfMat = Rcpp::as<MT>(T);
                ZP.linear_transformIt(TransorfMat.inverse());
            } else {
                TransorfMat = rotating < MT > (ZP, seed_rcpp);
            }
            Mat = extractMatPoly(ZP);
            break;
        }
        case 4: {
            throw Rcpp::exception("volesti does not support rotation for this representation currently.");
        }
    }

    TransorfMat.conservativeResize(n+1, n);
    TransorfMat.row(n) = VT::Ones(n);
    MT res(TransorfMat.rows(), Rcpp::as<MT>(Mat).rows()+n);
    res << Rcpp::as<MT>(Mat).transpose(), TransorfMat;
    return Rcpp::wrap(res);

}
