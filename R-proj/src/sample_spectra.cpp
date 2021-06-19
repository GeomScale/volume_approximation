// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 and 2019 program.

// Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "matrix_operations/EigenvaluesProblems.h"
#include "SDPAFormatManager.h"
#include "sampling/sampling.hpp"
#include "convex_bodies/spectrahedra_new/newSpectrahedron.h"



//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_spectra(Rcpp::Nullable<Rcpp::CharacterVector> file = R_NilValue,
                                  Rcpp::Nullable<unsigned int> N = R_NilValue,
                                  Rcpp::Nullable<unsigned int> walk_length = R_NilValue)
{
    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;

    std::list<Point> randPoints;

    SPECTRAHEDRON spectrahedron;
    Point objFunction;

    std::ifstream in;
    in.open(Rcpp::as<std::string>(file), std::ifstream::in);
    SdpaFormatManager<NT> sdpaFormatManager;
    sdpaFormatManager.loadSDPAFormatFile(in, spectrahedron, objFunction);

    std::vector<MT> mats = spectrahedron.getLMI().getMatrices();

    for (int i = 0; i < mats.size(); i++)
    {
        std::cout<<mats[i]<<"\n"<<std::endl;
    }

    std::cout<<objFunction.getCoefficients().transpose()<<std::endl;
    

    RNGType rng(spectrahedron.dimension());

    NT L = 0.2;
    unsigned int dim = spectrahedron.dimension();
    unsigned int walkL=1;
    unsigned int numpoints = 1000;
    unsigned int nburns = 0;
    if (walk_length.isNotNull()) {
        walkL = Rcpp::as<unsigned int>(walk_length);
    }
    if (N.isNotNull()) {
        numpoints = Rcpp::as<unsigned int>(N);
    }
    Point StartingPoint(dim);
    std::cout<<"dim = "<<dim<<std::endl;
    AcceleratedBilliardWalk WalkType(L);
    uniform_sampling(randPoints, spectrahedron, rng, WalkType, walkL, numpoints, StartingPoint, nburns);

    MT RetMat(dim, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++) {
        RetMat.col(jj) = (*rpit).getCoefficients();
    }

    return Rcpp::wrap(RetMat);
}

