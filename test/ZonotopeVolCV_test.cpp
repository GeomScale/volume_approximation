// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include "volume.h"
#include "polytope_generators.h"
#include "exact_vols.h"
#include <typeinfo>

template <typename NT>
NT factorial(NT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename NT>
void test_zono_volume(int n, int m, NT tolerance = 0.3)
{
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef Zonotope<Point> Zonotope;
    Zonotope ZP = gen_zonotope<Zonotope, RNGType>(n, m);

    // Setup the parameters
    int walk_len=1;
    int nexp=1, n_threads=1;
    NT e=0.1, err=0.0000000001;
    NT C=2.0,ratio,frac=0.1,delta=-1.0, round_value = 1.0;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    int N = 500 * ((int) C) + ((int) (n * n / 2));
    int W = 4*n*n+500;
    ratio = 1.0-1.0/(NT(n));
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);
    std::pair<NT,NT> res_round;

    vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,rng,
                          urdist,urdist1,-1.0,false,false,false,false,false,false,true,false);

    //Compute chebychev ball//
    std::pair<Point,NT> CheBall;
    //CheBall = ZP.ComputeInnerBall();

    // Estimate the volume
    std::cout << "--- Testing volume of Zonotope in dimension: " << n <<" and number of generators: "<< m << std::endl;
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    NT vol_exact = exact_zonotope_vol<NT>(ZP);
    //res_round = rounding_min_ellipsoid(ZP, CheBall, var);
    //round_value = round_value * res_round.first;
    NT vol = 0;
    unsigned int const num_of_exp = 10;
    for (unsigned int i=0; i<num_of_exp; i++)
    {
        CheBall = ZP.ComputeInnerBall();
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,rng,
                               urdist,urdist1,-1.0,false,false,false,false,false,false,true,false);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,CheBall.second,rng,C,frac,ratio,delta,false,
                                 false,false,false,false,false,false,true,false,false);
        vol += volume_gaussian_annealing(ZP, var1, var2, CheBall);
    }
    NT error = std::abs(((vol/num_of_exp)-vol_exact))/vol_exact;
    std::cout << "Computed volume (average) = " << vol/num_of_exp << std::endl;
    std::cout << "Expected volume = " << vol_exact << std::endl;
            CHECK(error < tolerance);
}


template <typename NT>
void call_test(int n, int m){
    test_zono_volume<NT>(n, m);
}


TEST_CASE("4_dimensional") {
    call_test<double>(4, 8);
    //call_test<float>(4,8);
    //call_test<long double>(4,8);
}

TEST_CASE("5_dimensional") {
    call_test<double>(5, 10);
    //call_test<float>(5, 10);
    //call_test<long double>(5, 10);
}
