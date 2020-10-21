#ifndef VOLESTIBINDINGS_H
#define VOLESTIBINDINGS_H

////// I think we could remove those...
//#include "Eigen/Eigen"
//#include "cartesian_kernel.h"

// from SOB volume - exactly the same for CG and CB methods
#include <fstream>
#include <iostream>
#include "random_walks.hpp"
#include "random.hpp"
#include "misc.h"
#include "known_polytope_generators.h"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"

//from generate_samples, some extra headers not already included
#include <chrono>
#include "sampling/sampling.hpp"

// for rounding
#include "preprocess/min_sampling_covering_ellipsoid_rounding.hpp"
#include "preprocess/svd_rounding.hpp"
#include "preprocess/svd_rounding_step.hpp"
#include "preprocess/max_inscribed_ellipsoid_rounding.hpp"
#include "preprocess/get_full_dimensional_polytope.hpp"

typedef double NT;
typedef Cartesian<NT>    Kernel;
typedef typename Kernel::Point    Point;
typedef HPolytope<Point> Hpolytope;
typedef typename Hpolytope::MT    MT;
typedef typename Hpolytope::VT    VT;
typedef BoostRandomNumberGenerator<boost::mt19937, double>    RNGType;


template <typename NT, typename MT, typename VT>
struct svd_rounding_parameters
{
public:

   svd_rounding_parameters() {}

    svd_rounding_parameters(unsigned int d)
        :   T(MT::Identity(d,d))
        ,   T_shift(VT::Zero(d))
        ,   num_rounding_steps(10*d)
        ,   fail(false)
        ,   converged(false)
        ,   last_round_under_p(false)
        ,   max_s(NT(0))
        ,   prev_max_s(std::numeric_limits<NT>::max())
        ,   round_it(1)
    {}

    MT T;
    VT T_shift;
    unsigned int num_rounding_steps;                                
    bool fail;
    bool converged;
    bool last_round_under_p;
    NT max_s;
    NT prev_max_s;
    unsigned int round_it;
};


// This is the HPolytopeCPP class; the main volesti class that is running the compute_volume(), rounding() and sampling() methods
class HPolytopeCPP{

   public:

      std::pair<Point,NT> CheBall;
      

      // regarding the rounding step
      typedef std::tuple<MT, VT, NT>    round_result;
      typedef svd_rounding_parameters<NT, MT, VT> svd_params;

      svd_params svd_parameters;

      // The class and its main specs
      HPolytopeCPP();
      HPolytopeCPP(double *A, double *b, int n_hyperplanes, int n_variables);

      Hpolytope HP;
      // Here we use the "~" destructor; this way we avoid a memory leak.
      ~HPolytopeCPP();

      // the compute_volume() function
      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed);

      // the generate_samples() function
      double generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, bool boundary, 
       bool cdhr, bool rdhr, bool gaussian, bool set_L, bool accelerated_billiard, bool billiard, bool ball_walk, double a, double L,  
       bool max_ball, double* inner_point, double radius,
       double* samples);

      // the rounding() function
      void rounding(char* rounding_method, double* new_A, double* new_b, double* T_matrix, double* shift, double &round_value,
       bool max_ball, double* inner_point, double radius);
      
      // the rounding_svd_step() function
      double rounding_svd_step(double* new_A, double* new_b, double* T_matrix, double* shift,
         double* inner_point, double radius);
      
};

// The preHPolytopeCPP class is responsible for the preprocess step of the polytope as well as for getting the full dimensional polytope
class lowDimHPolytopeCPP{

   public:

      MT A,Aeq;
      VT b,beq;

      // regarding getting full dimensional polytope
      typedef std::pair<Hpolytope, std::pair<MT, VT> > get_full_dim_pol_result;

      // The class and its main specs
      lowDimHPolytopeCPP();
      lowDimHPolytopeCPP(double *A, double *b, double *Aeq, double *beq, int n_rows_of_A, int n_cols_of_A, int n_row_of_Aeq, int n_cols_of_Aeq);
      Hpolytope low_HP;

      // Here we use the "~" destructor; this way we avoid a memory leak.
      ~lowDimHPolytopeCPP();

      // the get_full_dimensional_polytop() function
      int full_dimensiolal_polytope(double* N_extra_trans, double* shift, double* A_full_extra_trans, double* b_full);

};

#endif 

