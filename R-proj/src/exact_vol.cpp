// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/zpolytope.h"
#include "volume/exact_vols.h"

template <typename FT>
FT factorial(FT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

//' Compute the exact volume of (a) a zonotope (b) an arbitrary simplex in V-representation or (c) if the volume is known and declared by the input object.
//'
//' Given a zonotope (as an object of class Zonotope), this function computes the sum of the absolute values of the determinants of all the \eqn{d \times d} submatrices of the \eqn{m\times d} matrix \eqn{G} that contains row-wise the \eqn{m} \eqn{d}-dimensional segments that define the zonotope.
//' For an arbitrary simplex that is given in V-representation this function computes the absolute value of the determinant formed by the simplex's points assuming it is shifted to the origin.
//'
//' @param P A polytope
//'
//' @references \cite{E. Gover and N. Krikorian,
//' \dQuote{Determinants and the Volumes of Parallelotopes and Zonotopes,} \emph{Linear Algebra and its Applications, 433(1), 28 - 40,} 2010.}
//'
//' @return The exact volume of the input polytope, for zonotopes, simplices in V-representation and polytopes with known exact volume
//' @examples
//'
//' # compute the exact volume of a 5-dimensional zonotope defined by the Minkowski sum of 10 segments
//' Z = gen_rand_zonotope(2, 5)
//' vol = exact_vol(Z)
//'
//' \donttest{# compute the exact volume of a 2-d arbitrary simplex
//' V = matrix(c(2,3,-1,7,0,0),ncol = 2, nrow = 3, byrow = TRUE)
//' P = Vpolytope$new(V)
//' vol = exact_vol(P)
//' }
//'
//' # compute the exact volume the 10-dimensional cross polytope
//' P = gen_cross(10,'V')
//' vol = exact_vol(P)
//' @export
// [[Rcpp::export]]
double exact_vol(Rcpp::Nullable<Rcpp::Reference> P) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    if (NT(Rcpp::as<Rcpp::Reference>(P).field("volume")) > 0.0) {
        return NT(Rcpp::as<Rcpp::Reference>(P).field("volume"));
    }

    int type = Rcpp::as<Rcpp::Reference>(P).field("type"), dim;
    NT vol;

    if (type == 2) {

        dim = Rcpp::as<Rcpp::Reference>(P).field("dimension");

        if (Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).rows() ==
            Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).cols() + 1) {

            MT V = Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).transpose(), V2(dim,dim);
            VT v0 = V.col(dim);

            for (int i = 0; i < dim; ++i) {
                V2.col(i) = V.col(i) - v0;
            }
            vol = std::abs(V2.determinant()) / factorial(NT(dim));

        } else {
            throw Rcpp::exception("Volume unknown!");
        }

    } else if (type == 3) {

        typedef Zonotope<Point> zonotope;
        dim = Rcpp::as<Rcpp::Reference>(P).field("dimension");

        zonotope ZP(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")),
                VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")).rows()));
        vol = exact_zonotope_vol<NT>(ZP);

    } else {
        throw Rcpp::exception("Volume unknown!");
    }

    return vol;
}
