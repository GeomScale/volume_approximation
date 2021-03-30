// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef MMCS_HPP
#define MMCS_HPP

#include "diagnostics/ess_window_updater.hpp"

/**
   The class implements a single step of the Multiphase Monte Carlo Sampling algorithm
   given in,

   A. Chalkis, V. Fisikopoulos, E. Tsigaridas, H. Zafeiropoulos, Geometric algorithms for sampling the flux space of metabolic networks, SoCG 21.

 * @tparam Polytope convex polytope type
 * @tparam RandomNumberGenerator random number generator type
 * @tparam MT matrix type
 * @tparam Point cartensian point type
 * @tparam WalkTypePolicy random walk type
*/
template 
<
        typename Polytope,
        typename RandomNumberGenerator,
        typename MT,
        typename Point,
        typename WalkTypePolicy
>
void perform_mmcs_step(Polytope &P,
                       RandomNumberGenerator &rng,
                       const unsigned int &walk_length,
                       const unsigned int &target_ess,
                       unsigned int const& max_num_samples,
                       unsigned int const& window,
                       unsigned int &Neff_sampled,
                       unsigned int &total_samples,
                       unsigned int const& num_rounding_steps,
                       MT &TotalRandPoints,
                       bool &complete,
                       const Point &starting_point,
                       unsigned int const& nburns,
                       bool request_rounding,
                       WalkTypePolicy &WalkType)
{
    typedef typename Polytope::NT NT;
    typedef typename Polytope::VT VT;
    typedef typename WalkTypePolicy::template Walk
    <
        Polytope,
        RandomNumberGenerator
    > Walk;

    bool done = false;
    unsigned int points_to_sample = target_ess;
    int min_eff_samples;
    total_samples = 0;
    MT winPoints(P.dimension(), window);
    Point q(P.dimension());

    Point p = starting_point;

    if (request_rounding) 
    {
        TotalRandPoints.setZero(num_rounding_steps, P.dimension());
    } 
    else 
    {
        TotalRandPoints.setZero(max_num_samples, P.dimension());
    }
        
    Walk walk(P, p, rng, WalkType.param);
    ESSestimator<NT, VT, MT> estimator(window, P.dimension());

    walk.template parameters_burnin(P, p, 10 + int(std::log(NT(P.dimension()))), 10, rng);    

    while (!done) 
    {
        walk.template get_starting_point(P, p, q, 10, rng);
        for (int i = 0; i < window; i++)
        {
            walk.template apply(P, q, walk_length, rng);
            winPoints.col(i) = q.getCoefficients();
        }
        estimator.update_estimator(winPoints);
        total_samples += window;
        if (total_samples >= TotalRandPoints.rows()) 
        {
            if (total_samples > TotalRandPoints.rows()) 
            {
                TotalRandPoints.conservativeResize(total_samples, P.dimension());
            }
            if (request_rounding || total_samples >= max_num_samples)
            { 
                done = true;
            }
        }
        TotalRandPoints.block(total_samples - window, 0, window, P.dimension()) = winPoints.transpose();
        if (done || total_samples >= points_to_sample) 
        {
            estimator.estimate_effective_sample_size();                       
            min_eff_samples = int(estimator.get_effective_sample_size().minCoeff());
            if (done && min_eff_samples < target_ess) 
            {
                Neff_sampled = min_eff_samples;
                return;
            }
            if (min_eff_samples >= target_ess) 
            {
                complete = true;
                Neff_sampled = min_eff_samples;
                return;
            }
            if (min_eff_samples > 0) 
            {
                points_to_sample += (total_samples / min_eff_samples) * (target_ess - min_eff_samples) + 100;
            } 
            else 
            {
                points_to_sample = total_samples + 100;
            }
        }
    }
}

#endif

