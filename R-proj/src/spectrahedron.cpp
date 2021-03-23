// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "convex_bodies/spectrahedra/LMI.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "SDPAFormatManager.h"

//' Write a SDPA format file
//'
//' Outputs a spectrahedron (the matrices defining a linear matrix inequality) and a vector (the objective function)
//' to a SDPA format file.
//'
//' @param spectrahedron A spectrahedron in n dimensions; must be an object of class Spectrahedron
//' @param objectiveFunction A numerical vector of length n
//' @param outputFile Name of the output file
//'
//' @examples
//' \dontrun{
//' A0 = matrix(c(-1,0,0,0,-2,1,0,1,-2), nrow=3, ncol=3, byrow = TRUE)
//' A1 = matrix(c(-1,0,0,0,0,1,0,1,0), nrow=3, ncol=3, byrow = TRUE)
//' A2 = matrix(c(0,0,-1,0,0,0,-1,0,0), nrow=3, ncol=3, byrow = TRUE)
//' lmi = list(A0, A1, A2)
//' S = Spectrahedron$new(lmi);
//' objFunction = c(1,1)
//' writeSdpaFormatFile(S, objFunction, "output.txt")
//' }
//' @export
// [[Rcpp::export]]
void writeSdpaFormatFile(Rcpp::Nullable<Rcpp::Reference> spectrahedron = R_NilValue,
               Rcpp::Nullable<Rcpp::NumericVector> objectiveFunction = R_NilValue,
               Rcpp::Nullable<std::string> outputFile = R_NilValue) {

    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef LMI<NT, MT, VT> LMI;
    typedef Spectrahedron<NT, MT, VT> Spectrahedron;

    std::vector<MT> matrices = Rcpp::as<std::vector<MT> > (Rcpp::as<Rcpp::Reference> (spectrahedron).field("matrices"));
    LMI lmi(matrices);
    Spectrahedron _spectrahedron(lmi);
    Point c(Rcpp::as<VT> (objectiveFunction));

    std::filebuf fb;
    fb.open(Rcpp::as<std::string> (outputFile), std::ios::out);
    std::ostream os(&fb);

    SdpaFormatManager<NT> sdpaFormatManager;
    sdpaFormatManager.writeSDPAFormatFile(os, _spectrahedron, c);

    return;
}



//' Read a SDPA format file
//'
//' @param inputFile Name of the input file
//'
//' @return A list with two named items: an item "matrices" which is a list of the matrices and an vector "objFunction"
//'
//' @examples
//' path = system.file('extdata', package = 'volesti')
//' l = loadSdpaFormatFile(paste0(path,'/sdpa_n2m3.txt'))
//' @export
// [[Rcpp::export]]
Rcpp::List loadSdpaFormatFile(Rcpp::Nullable<std::string> inputFile = R_NilValue) {

    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef LMI<NT, MT, VT> LMI;
    typedef Spectrahedron<NT, MT, VT> Spectrahedron;

    Spectrahedron _spectrahedron;
    Point c;

    // open stream
    std::ifstream os;
    os.open(Rcpp::as<std::string> (inputFile),std::ifstream::in);

    // read file
    SdpaFormatManager<NT> sdpaFormatManager;
    sdpaFormatManager.loadSDPAFormatFile(os, _spectrahedron, c);

    std::vector<MT> matrices = _spectrahedron.getLMI().getMatrices();

    // return spectrahedron and objective function
    Rcpp::List _matrices;

    for (auto matrix : matrices)
        _matrices.push_back(Rcpp::wrap(matrix));

    Rcpp::List retList = Rcpp::List::create(Rcpp::Named("matrices") = _matrices , Rcpp::_["objFunction"] = Rcpp::wrap(c.getCoefficients()));
    return Rcpp::wrap(retList);
}
