// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef PARALLEL_MMCS_HPP
#define PARALLEL_MMCS_HPP


#include <iostream>
#include <omp.h>
#include <unistd.h>
#include "diagnostics/ess_window_updater.hpp"

/**
   The class implements a single step of the Parallel Multiphase Monte Carlo Sampling algorithm
   given in,

   A. Chalkis, V. Fisikopoulos, E. Tsigaridas, H. Zafeiropoulos, Geometric algorithms for sampling the flux space of metabolic networks, SoCG 21.

 * @tparam WalkTypePolicy random walk type
 * @tparam Polytope convex polytope type
 * @tparam RandomNumberGenerator random number generator type
 * @tparam MT matrix type
 * @tparam Point cartensian point type
 * @tparam NT number type
*/
template 
<
    typename WalkTypePolicy,
    typename Polytope,
    typename RandomNumberGenerator,
    typename MT,
    typename Point,
    typename NT
>
void perform_parallel_mmcs_step(Polytope &P,
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
                                unsigned int const& num_threads,
                                bool request_rounding,
                                NT L)
{
    typedef typename Polytope::VT VT; // vector type
    typedef typename WalkTypePolicy::template Walk
    <
        Polytope,
        RandomNumberGenerator
    > Walk;

    typedef typename WalkTypePolicy::template thread_parameters
    <
        NT,
        Point
    > _thread_parameters;

    omp_set_num_threads(num_threads);
    std::vector<int> num_starting_points_per_thread(num_threads, 0);
    std::vector<int> bound_on_num_points_per_thread(num_threads, 0);
    std::vector<int> num_generated_points_per_thread(num_threads, 0);
    unsigned int jj = 0, d = P.dimension(), m = P.num_of_hyperplanes();

    while (jj < nburns) 
    {
        for (unsigned int i = 0; i < num_threads; i++)
        {
            num_starting_points_per_thread[i]++;
            bound_on_num_points_per_thread[i] += window;
            jj++;
        }
    }

    std::vector<MT> winPoints_per_thread(num_threads, MT::Zero(d, window));
    std::vector<MT> TotalRandPoints_per_thread(num_threads);

    ESSestimator<NT, VT, MT> estimator(window, P.dimension());

    bool done = false, done_all = false;
    unsigned int points_to_sample = target_ess;
    int min_eff_samples;
    total_samples = 0;

    Point pp = starting_point;
    for (unsigned int i = 0; i < num_threads; i++)
    {
        TotalRandPoints_per_thread[i].setZero(bound_on_num_points_per_thread[i], d);
    }
    unsigned int upper_bound_on_total_num_of_samples;
    if (request_rounding) 
    {
        upper_bound_on_total_num_of_samples = num_rounding_steps;
    } 
    else 
    {
        upper_bound_on_total_num_of_samples = max_num_samples;
    }
    TotalRandPoints.resize(0, 0);
    Walk walk(P, L);
    
    _thread_parameters random_walk_parameters(d, m);
    walk.template parameters_burnin(P, pp, 10 + int(std::log(NT(d))), 10, rng, random_walk_parameters);
    Point const p = pp;

    #pragma omp parallel
    {
        int thread_index = omp_get_thread_num();
        _thread_parameters thread_random_walk_parameters(d, m);

        for (unsigned int it = 0; it < num_starting_points_per_thread[thread_index]; it++)    
        {
            if (done_all)
            {
                break;
            }
            walk.template get_starting_point(P, p, thread_random_walk_parameters, 10, rng);
            for (int i = 0; i < window; i++)
            {
                walk.template apply(P, thread_random_walk_parameters, walk_length, rng);
                winPoints_per_thread[thread_index].col(i) = thread_random_walk_parameters.p.getCoefficients();
            }
            estimator.update_estimator(winPoints_per_thread[thread_index]);
            num_generated_points_per_thread[thread_index] += window;
            total_samples += window;
            if (total_samples >= upper_bound_on_total_num_of_samples && thread_index == 0) 
            {
                done = true;
            }
            TotalRandPoints_per_thread[thread_index].block(num_generated_points_per_thread[thread_index] - window, 0, window, d) = winPoints_per_thread[thread_index].transpose();
            if ((done || total_samples >= points_to_sample) && thread_index == 0) 
            {
                estimator.estimate_effective_sample_size();                       
                min_eff_samples = int(estimator.get_effective_sample_size().minCoeff());
                if (done && min_eff_samples < target_ess) 
                {
                    Neff_sampled = min_eff_samples;
                    done_all = true;
                    continue;
                }
                if (min_eff_samples >= target_ess) 
                {
                    complete = true;
                    Neff_sampled = min_eff_samples;
                    done_all = true;
                    continue;
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

    estimator.estimate_effective_sample_size();                       
    min_eff_samples = int(estimator.get_effective_sample_size().minCoeff());
    Neff_sampled = min_eff_samples;
    if (min_eff_samples >= target_ess) 
    {
        complete = true;
    }

    unsigned int current_num_samples = 0;
    for (unsigned int i = 0; i < num_threads; i++)
    {
        TotalRandPoints.conservativeResize(TotalRandPoints.rows() + num_generated_points_per_thread[i], d);
        TotalRandPoints.block(current_num_samples, 0, num_generated_points_per_thread[i], d) = TotalRandPoints_per_thread[i].block(0, 0, num_generated_points_per_thread[i], d);
        current_num_samples += num_generated_points_per_thread[i];
        TotalRandPoints_per_thread[i].resize(0, 0);
    }
    
}

#endif

