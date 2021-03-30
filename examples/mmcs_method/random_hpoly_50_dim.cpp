// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

#include "Eigen/Eigen"

#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "sampling/mmcs.hpp"
#include "generators/h_polytopes_generator.h"
#include "diagnostics/multivariate_psrf.hpp"
#include "diagnostics/univariate_psrf.hpp"
#include "diagnostics/ess_window_updater.hpp"


template <typename NT>
void run_main() 
{
    typedef Cartesian<NT>    Kernel;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef boost::mt19937 PolyRNGType;
    typedef typename Kernel::Point    Point;
    typedef HPolytope <Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    int n = 50;

    RNGType rng(n);
    Hpolytope P = random_hpoly<Hpolytope, PolyRNGType>(n, 4*n, 127); // we fix the example polytope, seed = 127
    std::list<Point> randPoints;
    
    MT T = MT::Identity(n, n);
    VT T_shift = VT::Zero(n);

    unsigned int round_it = 1, num_rounding_steps = 20*n, 
                 walk_length = 1, num_its = 20, Neff = 4000, total_neff = 0, phase = 0,
                 window = 100, max_num_samples = 100 * n, total_samples, nburns, total_number_of_samples_in_P0 = 0;
    NT max_s, s_cutoff = 3.0, L;
    bool complete = false, request_rounding = true,
         rounding_completed = false;
    bool req_round_temp = request_rounding;

    std::pair<Point, NT> InnerBall;
    
    MT S;

    std::cout << "target effective sample size = " << Neff << "\n" << std::endl;
 
    while(true) 
    {
        phase++;
        if (request_rounding && rounding_completed) 
        {
            req_round_temp = false;
        }

        if (req_round_temp) 
        {
            nburns = num_rounding_steps / window + 1;
        } 
        else 
        {
            nburns = max_num_samples / window + 1;
        }

        InnerBall = P.ComputeInnerBall();
        L = NT(6) * std::sqrt(NT(n)) * InnerBall.second;
        AcceleratedBilliardWalk WalkType(L);

        unsigned int Neff_sampled;
        MT TotalRandPoints;
        perform_mmcs_step(P, rng, walk_length, Neff, max_num_samples, window, 
                          Neff_sampled, total_samples, num_rounding_steps, TotalRandPoints,
                          complete, InnerBall.first, 
                          nburns, req_round_temp, WalkType);

        Neff -= Neff_sampled;
        std::cout << "phase " << phase << ": number of correlated samples = " << total_samples << ", effective sample size = " << Neff_sampled;
        total_neff += Neff_sampled;
        Neff_sampled = 0;
        
        MT Samples = TotalRandPoints.transpose(); //do not copy TODO!
        for (int i = 0; i < total_samples; i++)
        {
            Samples.col(i) = T * Samples.col(i) + T_shift;
        }
        
        S.conservativeResize(P.dimension(), total_number_of_samples_in_P0 + total_samples);
        S.block(0, total_number_of_samples_in_P0, P.dimension(), total_samples) = Samples.block(0, 0, P.dimension(), total_samples);
        total_number_of_samples_in_P0 += total_samples;
        if (!complete) 
        {
            if (request_rounding && !rounding_completed) 
            {
                VT shift(n), s(n);
                MT V(n,n), S(n,n), round_mat;
                for (int i = 0; i < P.dimension(); ++i) 
                {
                    shift(i) = TotalRandPoints.col(i).mean();
                }

                for (int i = 0; i < total_samples; ++i) 
                {
                    TotalRandPoints.row(i) = TotalRandPoints.row(i) - shift.transpose();
                }

                Eigen::BDCSVD<MT> svd(TotalRandPoints, Eigen::ComputeFullV);
                s = svd.singularValues() / svd.singularValues().minCoeff();

                if (s.maxCoeff() >= 2.0) 
                {
                    for (int i = 0; i < s.size(); ++i) 
                    {
                        if (s(i) < 2.0) 
                        {
                            s(i) = 1.0;
                        }
                    }
                    V = svd.matrixV();
                } 
                else 
                {
                    s = VT::Ones(P.dimension());
                    V = MT::Identity(P.dimension(), P.dimension());
                }
                max_s = s.maxCoeff();
                S = s.asDiagonal();
                round_mat = V * S;

                round_it++;
                P.shift(shift);
                P.linear_transformIt(round_mat);
                T_shift += T * shift;
                T = T * round_mat;

                std::cout << ", ratio of the maximum singilar value over the minimum singular value = " << max_s << std::endl;

                if (max_s <= s_cutoff || round_it > num_its) 
                {
                    rounding_completed = true;
                }
            }
            else 
            {
                std::cout<<"\n";
            }
        } 
        else 
        {
            std::cout<<"\n\n";
            break;
        }
    } 

    std::cerr << "sum of effective sample sizes: " << total_neff << std::endl;
    std::cerr << "multivariate PSRF: " <<  multivariate_psrf<NT, VT>(S) << std::endl;
    std::cerr << "maximum marginal PSRF: " <<  univariate_psrf<NT, VT>(S).maxCoeff() << std::endl;
}

int main() {
  run_main<double>();
  return 0;
}
