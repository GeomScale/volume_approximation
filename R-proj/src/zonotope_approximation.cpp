// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "volume/volume_cooling_hpoly.hpp"

//' An internal Rccp function for the over-approximation of a zonotope
//'
//' @param Z A zonotope.
//' @param fit_ratio Optional. A boolean parameter to request the computation of the ratio of fitness.
//' @param settings Optional. A list that declares the values of the parameters of CB algorithm.
//' @param seed Optional. A fixed seed for the number generator.
//'
//' @keywords internal
//'
//' @return A List that contains a numerical matrix that describes the PCA approximation as a H-polytope and the ratio of fitness.
// [[Rcpp::export]]
Rcpp::List zono_approx (Rcpp::Reference Z,
                        Rcpp::Nullable<bool> fit_ratio = R_NilValue,
                        Rcpp::Nullable<Rcpp::List> settings = R_NilValue,
                        Rcpp::Nullable<double> seed = R_NilValue) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    typedef Zonotope <Point> zonotope;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    int n = Rcpp::as<int>(Z.field("dimension")), k = Rcpp::as<MT>(Z.field("G")).rows(), win_len = 200, walkL = 1;

    RNGType rng(n);
    if (seed.isNotNull()) {
        unsigned seed_rcpp = Rcpp::as<double>(seed);
        rng.set_seed(seed_rcpp);
    }

    NT e = 0.1, ratio = std::numeric_limits<double>::signaling_NaN();
    bool hpoly = false;

    MT X(n, 2 * k);
    X << Rcpp::as<MT>(Z.field("G")).transpose(), -Rcpp::as<MT>(Z.field("G")).transpose();
    Eigen::JacobiSVD <MT> svd(X * X.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
    MT G(k, 2 * n);
    G << Rcpp::as<MT>(Z.field("G")) * svd.matrixU(), Rcpp::as<MT>(Z.field("G")) * svd.matrixU();
    VT Gred_ii = G.transpose().cwiseAbs().rowwise().sum();
    MT A(n, 2 * n);
    A << -MT::Identity(n, n), MT::Identity(n, n);
    MT Mat(2 * n, n + 1);
    Mat << Gred_ii, A.transpose() * svd.matrixU().transpose();

    //Hpolytope HP(n, A.transpose() * svd.matrixU().transpose(), Gred_ii);

    if (fit_ratio.isNotNull() && Rcpp::as<bool>(fit_ratio)) {
        NT vol_red = std::abs(svd.matrixU().determinant());
        for (int i = 0; i < n; ++i) {
            vol_red *= 2.0 * Gred_ii(i);
        }

        walkL = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("walk_length")) ? 1 : Rcpp::as<int>(
                Rcpp::as<Rcpp::List>(settings)["walk_length"]);
        e = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("error")) ? 0.1 : Rcpp::as<NT>(
                Rcpp::as<Rcpp::List>(settings)["error"]);
        win_len = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("win_len")) ? 200 : Rcpp::as<int>(
                Rcpp::as<Rcpp::List>(settings)["win_len"]);

        zonotope ZP(n, Rcpp::as<MT>(Z.field("G")), VT::Ones(Rcpp::as<MT>(Z.field("G")).rows()));

        if (Rcpp::as<Rcpp::List>(settings).containsElementNamed("hpoly")) {
            hpoly = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(settings)["hpoly"]);
        } else if (ZP.num_of_generators() / ZP.dimension() < 5 ) {
            hpoly = true;
        } else {
            hpoly = false;
        }

        NT vol;
        if (!hpoly) {
            vol = volume_cooling_balls<BilliardWalk>(ZP, rng, e, walkL, win_len);
        } else {
            vol = volume_cooling_hpoly<BilliardWalk, Hpolytope>(ZP, rng, e, walkL, win_len);
        }
        ratio = std::pow(vol_red / vol, 1.0 / NT(n));
    }

    return Rcpp::List::create(Rcpp::Named("Mat") = Rcpp::wrap(Mat), Rcpp::Named("fit_ratio") = ratio);

}
