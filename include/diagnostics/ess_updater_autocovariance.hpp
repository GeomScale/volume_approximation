// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef DIAGNOSTICS_ESS_UPDATER_AUTOCOVARIANCE_HPP
#define DIAGNOSTICS_ESS_UPDATER_AUTOCOVARIANCE_HPP

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <unsupported/Eigen/FFT>
#include <complex>
#include <cmath>


/**
    Compute autocovariance estimates for every lag for the input sequence
    using the Geyer's stable estimator given in 
    Charles J. Geyer, Practical Markov Chain Monte Carlo, Statistical Science 1992.

    * @tparam NT number type
    * @tparam VT vector type
    * @param samples the sequence of correlated samples
    * @param auto_cov the autocovariance estimates
*/
template <typename NT, typename VT>
void compute_autocovariance(VT const& samples, VT &auto_cov) 
{
    typedef Eigen::FFT<NT> EigenFFT;
    typedef Eigen::Matrix<std::complex<NT>, Eigen::Dynamic, 1> CVT;
    EigenFFT fft;

    unsigned int N = samples.size();
    NT samples_mean = samples.mean();
    auto_cov.setZero(N);

    // compute normalized samples
    VT normalized_samples(2 * N);
    normalized_samples.setZero();
    normalized_samples.head(N) = samples.array() - samples_mean;

    NT variance = (normalized_samples.cwiseProduct(normalized_samples)).sum();
    variance *= (1.0 / N);
    variance += NT(1e-16)*(samples_mean*samples_mean);
    normalized_samples.head(N) = normalized_samples.head(N).array() / sqrt(variance);

    // Perform FFT on 2N points
    CVT frequency(2 * N);
    fft.fwd(frequency, normalized_samples);

    // Invert fft to get autocorrelation function
    CVT auto_cov_tmp(2 * N);
    frequency = frequency.cwiseAbs2();
    fft.inv(auto_cov_tmp, frequency);

    auto_cov = auto_cov_tmp.head(N).real().array() / N;

    boost::accumulators::accumulator_set<NT, boost::accumulators::stats<boost::accumulators::tag::variance>> accumulator;
    for (int i = 0; i < samples.size(); ++i) 
    {
        accumulator(samples.coeff(i));
    }
  
    auto_cov = auto_cov.array() * boost::accumulators::variance(accumulator);
}

#endif
