// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 and 2019 program.


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "volume/volume_sequence_of_balls.hpp"
#include "sampling/simplex.hpp"


//' Sample perfect uniformly distributed points from well known convex bodies: (a) the unit simplex, (b) the canonical simplex, (c) the boundary of a hypersphere or (d) the interior of a hypersphere.
//'
//' The \eqn{d}-dimensional unit simplex is the set of points \eqn{\vec{x}\in \R^d}, s.t.: \eqn{\sum_i x_i\leq 1}, \eqn{x_i\geq 0}. The \eqn{d}-dimensional canonical simplex is the set of points \eqn{\vec{x}\in \R^d}, s.t.: \eqn{\sum_i x_i = 1}, \eqn{x_i\geq 0}.
//'
//' @param body A list to request exact uniform sampling from special well known convex bodies through the following input parameters:
//' \itemize{
//' \item{\code{type} }{ A string that declares the type of the body for the exact sampling: a) \code{'unit_simplex'} for the unit simplex, b) \code{'canonical_simplex'} for the canonical simplex, c) \code{'hypersphere'} for the boundary of a hypersphere centered at the origin, d) \code{'ball'} for the interior of a hypersphere centered at the origin.}
//' \item{\code{dimension} }{ An integer that declares the dimension when exact sampling is enabled for a simplex or a hypersphere.}
//' \item{\code{radius} }{ The radius of the \eqn{d}-dimensional hypersphere. The default value is \eqn{1}.}
//' }
//' @param n The number of points that the function is going to sample.
//' @param seed Optional. A fixed seed for the number generator.
//'
//' @references \cite{R.Y. Rubinstein and B. Melamed,
//' \dQuote{Modern simulation and modeling} \emph{ Wiley Series in Probability and Statistics,} 1998.}
//' @references \cite{A Smith, Noah and W Tromble, Roy,
//' \dQuote{Sampling Uniformly from the Unit Simplex,} \emph{ Center for Language and Speech Processing Johns Hopkins University,} 2004.}
//'
//' @return A \eqn{d\times n} matrix that contains, column-wise, the sampled points from the convex polytope P.
//' @examples
//' # 100 uniform points from the 2-d unit ball
//' points = direct_sampling(n = 100, body = list("type" = "ball", "dimension" = 2))
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix direct_sampling(Rcpp::Nullable<Rcpp::List> body,
                                    Rcpp::Nullable<unsigned int> n,
                                    Rcpp::Nullable<double> seed = R_NilValue) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef boost::mt19937 RNGType2;
    typedef typename Kernel::Point Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;

    int dim, numpoints;
    NT radius = 1.0;
    std::list<Point> randPoints;

    if (!Rcpp::as<Rcpp::List>(body).containsElementNamed("dimension")) {
        throw Rcpp::exception("Dimension has to be given as input!");
    }
    dim = Rcpp::as<int>(Rcpp::as<Rcpp::List>(body)["dimension"]);
    if (dim <=1) throw Rcpp::exception("Dimension has to be larger than 1!");
    RNGType rng(dim);

    if (seed.isNotNull()) {
        unsigned seed2 = Rcpp::as<double>(seed);
        rng.set_seed(seed2);
    }
    double seed3 = (!seed.isNotNull()) ? std::numeric_limits<double>::signaling_NaN() : Rcpp::as<double>(seed);
    //RNGType rng2(5);

    numpoints = (!n.isNotNull()) ? 100 : Rcpp::as<unsigned int>(n);

    numpoints = Rcpp::as<unsigned int>(n);
    if (numpoints <= 0) throw Rcpp::exception("The number of samples has to be a positice integer!");


    if (Rcpp::as<Rcpp::List>(body).containsElementNamed("radius")) {

        radius = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(body)["radius"]);
        if (radius <= NT(0)) throw Rcpp::exception("Radius has to be a positive number!");

    }
    if (!Rcpp::as<Rcpp::List>(body).containsElementNamed("type")) {

        throw Rcpp::exception("The kind of body has to be given as input!");

    }
    if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(body)["type"]).compare(std::string("hypersphere"))==0) {

        for (unsigned int k = 0; k < numpoints; ++k) {
            randPoints.push_back(GetPointOnDsphere<Point>::apply(dim, radius, rng));
        }

    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(body)["type"]).compare(std::string("ball"))==0) {

        for (unsigned int k = 0; k < numpoints; ++k) {
            randPoints.push_back(GetPointInDsphere<Point>::apply(dim, radius, rng));
        }

    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(body)["type"]).compare(std::string("unit_simplex"))==0) {

        Sam_Unit<NT, RNGType2 >(dim, numpoints, randPoints, seed3);

    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(body)["type"]).compare(std::string("canonical_simplex"))==0) {

        Sam_Canon_Unit<NT, RNGType2 >(dim, numpoints, randPoints, seed3);

    } else {

        throw Rcpp::exception("Wrong input!");

    }

    MT RetMat(dim, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++) {
        RetMat.col(jj) = rpit->getCoefficients();
    }
    return Rcpp::wrap(RetMat);
}
