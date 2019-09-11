// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include <fstream>
#include "volume.h"
#include "polytope_generators.h"
#include <string>
#include <typeinfo>

template <typename NT>
NT factorial(NT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename NT, class RNGType, class Polytope>
void test_CV_volume(Polytope &HP, NT expected, NT tolerance=0.3)
{

    typedef typename Polytope::PolytopePoint Point;

    // Setup the parameters
    int n = HP.dimension();
    int walk_len=1;
    int nexp=1, n_threads=1;
    NT e=0.1, err=0.0000000001;
    NT C=2.0,ratio,frac=0.1,delta=-1.0;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    int N = 500 * ((int) C) + ((int) (n * n / 2));
    int W = 4*n*n+500;
    ratio = 1.0-1.0/(NT(n));
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    //compute the chebychev ball
    std::pair<Point,NT> CheBall;


    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    NT vol = 0;
    unsigned int const num_of_exp = 20;
    for (unsigned int i=0; i<num_of_exp; i++)
    {
        CheBall = HP.ComputeInnerBall();
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,true,false);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,CheBall.second,rng,C,frac,ratio,delta,false,
                    false,false,false,false,false,false,true,false,false);
        vol += volume_gaussian_annealing(HP, var1, var2, CheBall);
    }
    NT error = std::abs(((vol/num_of_exp)-expected))/expected;
    std::cout << "Computed volume (average) = " << vol/num_of_exp << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
            CHECK(error < tolerance);
}


template <typename NT>
void call_test_cube(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-cube10" << std::endl;
    P = gen_cube<Hpolytope>(10, false);
    test_CV_volume<NT, RNGType>(P, 1024.0);

    std::cout << "--- Testing volume of H-cube20" << std::endl;
    P = gen_cube<Hpolytope>(20, false);
    test_CV_volume<NT, RNGType>(P, 1048576.0);

    std::cout << "--- Testing volume of H-cube30" << std::endl;
    P = gen_cube<Hpolytope>(30, false);
    test_CV_volume<NT, RNGType>(P, 1073742000.0, 0.3);

}

template <typename NT>
void call_test_cross(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;

    std::cout << "--- Testing volume of H-cross10" << std::endl;
    Hpolytope P = gen_cross<Hpolytope>(10, false);
    test_CV_volume<NT, RNGType>(P, 0.0002821869);

}

template <typename NT>
void call_test_birk() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-birk3" << std::endl;
    std::ifstream inp;
    std::vector<std::vector<NT> > Pin;
    inp.open("../R-proj/inst/extdata/birk3.ine",std::ifstream::in);
    read_pointset(inp,Pin);
    P.init(Pin);
    test_CV_volume<NT, RNGType>(P, 0.125);

    //test_CV_volume<NT>("../data/birk4.ine", 0.000970018);
    //test_CV_volume<NT>("../data/birk5.ine", 0.000000225);

    //std::cout << "--- Testing volume of H-birk6" << std::endl;
    //std::ifstream inp4;
    //std::vector<std::vector<NT> > Pin4;
    //inp4.open("../R-proj/inst/extdata/birk6.ine",std::ifstream::in);
    //read_pointset(inp4,Pin4);
    //P.init(Pin4);
    //test_CV_volume<NT, RNGType>(P, 0.0000000000009455459196, 0.5);

}

template <typename NT>
void call_test_prod_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-prod_simplex5" << std::endl;
    P = gen_prod_simplex<Hpolytope>(5);
    test_CV_volume<NT, RNGType>(P, std::pow(1.0 / factorial(5.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex10" << std::endl;
    P = gen_prod_simplex<Hpolytope>(10);
    test_CV_volume<NT, RNGType>(P, std::pow(1.0 / factorial(10.0), 2));

    //std::cout << "--- Testing volume of H-prod_simplex15" << std::endl;
    //P = gen_prod_simplex<Hpolytope>(15);
    //test_CV_volume<NT, RNGType>(P, std::pow(1.0 / factorial(15.0), 2), 0.3);

    //std::cout << "--- Testing volume of H-prod_simplex20" << std::endl;
    //P = gen_prod_simplex<Hpolytope>(20);
    //test_CV_volume<NT, RNGType>(P, std::pow(1.0 / factorial(20.0), 2), 0.3);

}

template <typename NT>
void call_test_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-simplex10" << std::endl;
    P = gen_simplex<Hpolytope>(10, false);
    test_CV_volume<NT, RNGType>(P, 1.0 / factorial(10.0), 0.2);

    std::cout << "--- Testing volume of H-simplex20" << std::endl;
    P = gen_simplex<Hpolytope>(20, false);
    test_CV_volume<NT, RNGType>(P, 1.0 / factorial(20.0), 0.2);

    std::cout << "--- Testing volume of H-simplex30" << std::endl;
    P = gen_simplex<Hpolytope>(30, false);
    test_CV_volume<NT, RNGType>(P, 1.0 / factorial(30.0), 0.3);

    //std::cout << "--- Testing volume of H-simplex40" << std::endl;
    //P = gen_simplex<Hpolytope>(40, false);
    //test_CV_volume<NT, RNGType>(P, 1.0 / factorial(40.0), 0.3);

    //std::cout << "--- Testing volume of H-simplex50" << std::endl;
    //P = gen_simplex<Hpolytope>(50, false);
    //test_CV_volume<NT, RNGType>(P, 1.0 / factorial(50.0), 0.3);

}


TEST_CASE("cube") {
    call_test_cube<double>();
    //call_test_cube<float>();
    //call_test_cube<long double>();
}

TEST_CASE("cross") {
    call_test_cross<double>();
    //call_test_cross<float>();
    //call_test_cross<long double>();
}

TEST_CASE("birk") {
    call_test_birk<double>();
    //call_test_birk<float>();
    //call_test_birk<long double>();
}

TEST_CASE("prod_simplex") {
    call_test_prod_simplex<double>();
    //call_test_prod_simplex<float>();
    //call_test_prod_simplex<long double>();
}

TEST_CASE("simplex") {
    call_test_simplex<double>();
    //call_test_simplex<float>();
    //call_test_simplex<long double>();
}
