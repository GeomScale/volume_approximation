// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "diagnostics/raftery.hpp"

//' Raftery and Lewis MCMC diagnostic
//'
//' @param samples A matrix that contans column-wise the sampled points from a geometric random walk.
//' @param q Optional. The quantile of the quantity of interest. The default value is 0.025.
//' @param r Optional. The level of precision desired. The default value is 0.01.
//' @param s Optional. The probability associated with r. The default value is 0.95.
//'
//' @references \cite{Raftery, A. E. and Lewis, S. M.,
//' \dQuote{How many iterations in the Gibbs sampler?,} \emph{Bayesian Statistics 4. Proceedings of the Fourth Valencia International Meeting,} 1992.}
//'
//' @return (i) The number of draws required for burn-in, (ii) the skip parameter for 1st-order Markov chain, (iii) the skip parameter sufficient to get independence chain, (iv) the number of draws required to achieve r precision, (v) the number of draws if the chain is white noise, (vi) the I-statistic from Raftery and Lewis (1992).
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix raftery(Rcpp::NumericMatrix samples,
                            Rcpp::Nullable<double> q = R_NilValue,
                            Rcpp::Nullable<double> r = R_NilValue,
                            Rcpp::Nullable<double> s = R_NilValue)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    NT _q = (!q.isNotNull()) ? 0.025 : Rcpp::as<NT>(q);
    NT _r = (!r.isNotNull()) ? 0.01 : Rcpp::as<NT>(r);
    NT _s = (!s.isNotNull()) ? 0.95 : Rcpp::as<NT>(s);

    MT results = perform_raftery<VT>(Rcpp::as<MT>(samples), _q, _r, _s);

    return Rcpp::wrap(results);
}
