// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20014-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "diagnostics/geweke.hpp"


//' Geweke's MCMC diagnostic
//'
//' @param samples A matrix that contans column-wise the sampled points from a geometric random walk.
//' @param frac_first Optional. The portion of the first in order points in matrix samples.
//' @param frac_last Optional. The portion of the last in order points in matrix samples.
//'
//' @references \cite{Geweke, J.,
//' \dQuote{Evaluating the accuracy of sampling-based approaches to the calculation of posterior moments,} \emph{ In Bayesian Statistics 4. Proceedings of the Fourth Valencia International Meeting,} 1992.}
//'
//' @return A boolean to denote if the result of Geweke diagnostic: (i)  false if the null hypothesis is rejected, (ii) true if the null hypothesis is not rejected.
//'
//' @export
// [[Rcpp::export]]
bool geweke(Rcpp::NumericMatrix samples, 
            Rcpp::Nullable<double> frac_first = R_NilValue, 
            Rcpp::Nullable<double> frac_last = R_NilValue)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    NT frac_1, frac_2;

    frac_1 = (!frac_first.isNotNull()) ? NT(0.1) : Rcpp::as<NT>(frac_first);
    frac_2 = (!frac_last.isNotNull()) ? NT(0.5) : Rcpp::as<NT>(frac_last);

    return perform_geweke<VT>(Rcpp::as<MT>(samples), frac_1, frac_2);
}
