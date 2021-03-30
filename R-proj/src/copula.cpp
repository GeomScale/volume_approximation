// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "convex_bodies/ellipsoids.h"
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include "sampling/simplex.hpp"
#include "volume/copulas.h"

//' Construct a copula using uniform sampling from the unit simplex
//'
//' Given two families of parallel hyperplanes or a family of parallel hyperplanes and a family of concentric ellispoids centered at the origin intersecting the canonical simplex, this function uniformly samples from the canonical simplex and construct an approximation of the bivariate probability distribution, called copula (see \url{https://en.wikipedia.org/wiki/Copula_(probability_theory)}).
//' At least two families of hyperplanes or one family of hyperplanes and one family of ellipsoids have to be given as input.
//'
//' @param r1 The \eqn{d}-dimensional normal vector of the first family of parallel hyperplanes.
//' @param r2 Optional. The \eqn{d}-dimensional normal vector of the second family of parallel hyperplanes.
//' @param sigma Optional. The \eqn{d\times d} symmetric positive semidefine matrix that describes the family of concentric ellipsoids centered at the origin.
//' @param m The number of the slices for the copula. The default value is 100.
//' @param n The number of points to sample. The default value is \eqn{5\cdot 10^5}.
//' @param seed Optional. A fixed seed for the number generator.
//'
//' @references \cite{L. Cales, A. Chalkis, I.Z. Emiris, V. Fisikopoulos,
//' \dQuote{Practical volume computation of structured convex bodies, and an application to modeling portfolio dependencies and financial crises,} \emph{Proc. of Symposium on Computational Geometry, Budapest, Hungary,} 2018.}
//'
//' @return A \eqn{m\times m} numerical matrix that corresponds to a copula.
//' @examples
//' # compute a copula for two random families of parallel hyperplanes
//' h1 = runif(n = 10, min = 1, max = 1000)
//' h1 = h1 / 1000
//' h2=runif(n = 10, min = 1, max = 1000)
//' h2 = h2 / 1000
//' cop = copula(r1 = h1, r2 = h2, m = 10, n = 100000)
//'
//' # compute a copula for a family of parallel hyperplanes and a family of conentric ellipsoids
//' h = runif(n = 10, min = 1, max = 1000)
//' h = h / 1000
//' E = replicate(10, rnorm(20))
//' E = cov(E)
//' cop = copula(r1 = h, sigma = E, m = 10, n = 100000)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix copula (Rcpp::Nullable<Rcpp::NumericVector> r1,
                            Rcpp::Nullable<Rcpp::NumericVector> r2 = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericMatrix> sigma = R_NilValue,
                            Rcpp::Nullable<unsigned int> m = R_NilValue,
                            Rcpp::Nullable<unsigned int> n = R_NilValue,
                            Rcpp::Nullable<double> seed = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef copula_ellipsoid<Point, MT, VT> CopEll;
    unsigned int num_slices = 100, numpoints = 500000;

    if (m.isNotNull()) {
        num_slices = Rcpp::as<unsigned int>(m);
    }

    if (n.isNotNull()) {
        numpoints = Rcpp::as<unsigned int>(n);
    }

    double seed3 = (!seed.isNotNull()) ? std::numeric_limits<double>::signaling_NaN() : Rcpp::as<double>(seed);

    Rcpp::NumericMatrix copula(num_slices, num_slices);
    std::vector<std::vector<NT> > StdCopula;
    unsigned int dim = Rcpp::as<std::vector<NT> >(r1).size(), i, j;

    std::vector<NT> hyp1 = Rcpp::as<std::vector<NT> >(r1);
    if (r2.isNotNull()) {

        std::vector <NT> hyp2 = Rcpp::as < std::vector < NT > > (r2);
        StdCopula = twoParHypFam<Point, RNGType>(dim, numpoints, num_slices, hyp1, hyp2, seed3);

    } else if (sigma.isNotNull()) {

        std::vector<std::vector<NT> > Gin(dim, std::vector<NT>(dim));
        MT EE = Rcpp::as<MT>(sigma);
        for (int i=0; i<dim; i++) {
            for (int j = 0; j < dim; j++) {
                Gin[i][j] = EE(i, j);
            }
        }
        CopEll Ell(Gin);
        StdCopula = hypfam_ellfam<Point, RNGType >(dim, numpoints, num_slices, hyp1, Ell, seed3);
    } else {

        throw Rcpp::exception("You have to give as input either two families of hyperplanes or one family of hyperplanes and one family of ellipsoids!");

    }

    for(int i=0; i<num_slices; i++) {
        for(int j=0; j<num_slices; j++){
            copula(i,j) = StdCopula[i][j];
        }
    }

    return copula;
}
