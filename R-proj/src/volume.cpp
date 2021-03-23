// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 and 2019 program.


#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "random_walks/random_walks.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "volume/volume_cooling_hpoly.hpp"
#include "preprocess/max_inscribed_ellipsoid_rounding.hpp"
#include "preprocess/svd_rounding.hpp"

enum random_walks {ball_walk, rdhr, cdhr, billiard, accelarated_billiard};
enum volume_algorithms {CB, CG, SOB};
enum rounding_type {none, min_ellipsoid, max_ellipsoid, svd};

template <typename Polytope,  typename RNGType,  typename NT>
double generic_volume(Polytope& P, RNGType &rng, unsigned int walk_length, NT e,
                      volume_algorithms const& algo, unsigned int win_len,
                      rounding_type const& rounding, random_walks const& walk)
{
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;
    typedef typename Polytope::PointType Point;

    typedef HPolytope<Point> Hpolytope;

    NT round_val = 1.0;
    unsigned int n = P.dimension();
    std::pair<Point, NT> InnerBall;

    if (rounding != none){
         InnerBall = P.ComputeInnerBall();
         if (InnerBall.second < 0.0) throw Rcpp::exception("Unable to compute a feasible point.");
    }

    switch (rounding)
    {
    case min_ellipsoid:
        switch (walk)
        {
        case cdhr:
            round_val = std::get<2>(min_sampling_covering_ellipsoid_rounding<CDHRWalk, MT, VT>(P, InnerBall, 10 + 10 * n, rng));
            break;
        case accelarated_billiard:
            round_val = std::get<2>(min_sampling_covering_ellipsoid_rounding<AcceleratedBilliardWalk, MT, VT>(P, InnerBall, 2, rng));
            break;
        default:
            round_val = std::get<2>(min_sampling_covering_ellipsoid_rounding<BilliardWalk, MT, VT>(P, InnerBall, 2, rng));
            break;
        }
        break;
    case svd:
        switch (walk)
        {
        case cdhr:
            round_val = std::get<2>(svd_rounding<CDHRWalk, MT, VT>(P, InnerBall, 10 + 10 * n, rng));
            break;
        case accelarated_billiard:
            round_val = std::get<2>(svd_rounding<AcceleratedBilliardWalk, MT, VT>(P, InnerBall, 2, rng));
            break;
        default:
            round_val = std::get<2>(svd_rounding<BilliardWalk, MT, VT>(P, InnerBall, 2, rng));
            break;
        }
        break;
    case max_ellipsoid:
        round_val = std::get<2>(max_inscribed_ellipsoid_rounding<MT, VT, NT>(P, InnerBall.first));
        break;
    default:
        break;
    }

    NT vol;
    switch (algo)
    {
    case CG:
        switch (walk)
        {
        case cdhr:
            vol = volume_cooling_gaussians<GaussianCDHRWalk>(P, rng, e, walk_length);
            break;
        case rdhr:
            vol = volume_cooling_gaussians<GaussianRDHRWalk>(P, rng, e, walk_length);
            break;
        case ball_walk:
            vol = volume_cooling_gaussians<GaussianBallWalk>(P, rng, e, walk_length);
            break;
        default:
            throw Rcpp::exception("This random walk can not be used by CG algorithm!");
            break;
        }
        break;
    case CB:
        switch (walk)
        {
        case cdhr:
            vol = volume_cooling_balls<CDHRWalk>(P, rng, e, walk_length, win_len);
            break;
        case rdhr:
            vol = volume_cooling_balls<RDHRWalk>(P, rng, e, walk_length, win_len);
            break;
        case ball_walk:
            vol = volume_cooling_balls<BallWalk>(P, rng, e, walk_length, win_len);
            break;
        case billiard:
            vol = volume_cooling_balls<BilliardWalk>(P, rng, e, walk_length, win_len);
            break;
        case accelarated_billiard:
            vol = volume_cooling_balls<AcceleratedBilliardWalk>(P, rng, e, walk_length, win_len);
            break;
        default:
            throw Rcpp::exception("This random walk can not be used by CB algorithm!");
            break;
        }
        break;
    case SOB:
        switch (walk)
        {
        case cdhr:
            vol = volume_sequence_of_balls<CDHRWalk>(P, rng, e, walk_length);
            break;
        case rdhr:
            vol = volume_sequence_of_balls<RDHRWalk>(P, rng, e, walk_length);
            break;
        case ball_walk:
            vol = volume_sequence_of_balls<BallWalk>(P, rng, e, walk_length);
            break;
        case billiard:
            vol = volume_sequence_of_balls<BilliardWalk>(P, rng, e, walk_length);
            break;
        case accelarated_billiard:
            vol = volume_sequence_of_balls<AcceleratedBilliardWalk>(P, rng, e, walk_length);
            break;
        default:
            throw Rcpp::exception("This random walk can not be used by CB algorithm!");
            break;
        }
        break;
    default:
        throw Rcpp::exception("Unknown algorithm!");
        break;
    }
    if (vol < 0.0) throw Rcpp::exception("volesti failed to terminate.");
    vol *= round_val;
    return vol;
}

//' The main function for volume approximation of a convex Polytope (H-polytope, V-polytope, zonotope or intersection of two V-polytopes)
//'
//' For the volume approximation can be used three algorithms. Either CoolingBodies (CB) or SequenceOfBalls (SOB) or CoolingGaussian (CG). An H-polytope with \eqn{m} facets is described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{P=\{x\ |\  Ax\leq b\} }. A V-polytope is defined as the convex hull of \eqn{m} \eqn{d}-dimensional points which correspond to the vertices of P. A zonotope is desrcibed by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
//'
//' @param P A convex polytope. It is an object from class a) Hpolytope or b) Vpolytope or c) Zonotope or d) VpolytopeIntersection.
//' @param settings Optional. A list that declares which algorithm, random walk and values of parameters to use, as follows:
//' \itemize{
//' \item{\code{algorithm} }{ A string to set the algorithm to use: a) \code{'CB'} for CB algorithm, b) \code{'SoB'} for SOB algorithm or b) \code{'CG'} for CG algorithm. The defalut algorithm is \code{'CB'}.}
//' \item{\code{error} }{ A numeric value to set the upper bound for the approximation error. The default value is \eqn{1} for SOB algorithm and \eqn{0.1} otherwise.}
//' \item{\code{random_walk} }{ A string that declares the random walk method: a) \code{'CDHR'} for Coordinate Directions Hit-and-Run, b) \code{'RDHR'} for Random Directions Hit-and-Run, c) \code{'BaW'} for Ball Walk, or \code{'BiW'} for Billiard walk. For CB algorithm the default walk is \code{'BiW'}. For CG and SOB algorithms the default walk is \code{'CDHR'} for H-polytopes and \code{'RDHR'} for the other representations.}
//' \item{\code{walk_length} }{ An integer to set the number of the steps for the random walk. The default value is \eqn{\lfloor 10 + d/10\rfloor} for \code{'SOB'} and \eqn{1} otherwise.}
//' \item{\code{win_len} }{ The length of the sliding window for CB or CG algorithm. The default value is \eqn{250} for CB with BiW and \eqn{400+3d^2} for CB and any other random walk and \eqn{500+4d^2} for CG.}
//' \item{\code{hpoly} }{ A boolean parameter to use H-polytopes in MMC of CB algorithm when the input polytope is a zonotope. The default value is \code{TRUE} when the order of the zonotope is \eqn{<5}, otherwise it is \code{FALSE}.}
//' }
//' @param rounding Optional. A string parameter to request a rounding method to be applied in the input polytope before volume computation: a) \code{'min_ellipsoid'}, b) \code{'svd'}, c) \code{'max_ellipsoid'} and d) \code{'none'} for no rounding.
//' @param seed Optional. A fixed seed for the number generator.
//'
//' @references \cite{I.Z.Emiris and V. Fisikopoulos,
//' \dQuote{Practical polytope volume approximation,} \emph{ACM Trans. Math. Soft.,} 2018.},
//' @references \cite{A. Chalkis and I.Z.Emiris and V. Fisikopoulos,
//' \dQuote{Practical Volume Estimation by a New Annealing Schedule for Cooling Convex Bodies,} \emph{CoRR, abs/1905.05494,} 2019.},
//' @references \cite{B. Cousins and S. Vempala, \dQuote{A practical volume algorithm,} \emph{Springer-Verlag Berlin Heidelberg and The Mathematical Programming Society,} 2015.}
//'
//'
//' @return The approximation of the volume of a convex polytope.
//' @examples
//'
//' # calling SOB algorithm for a H-polytope (5d unit simplex)
//' HP = gen_cube(5,'H')
//' vol = volume(HP)
//'
//' # calling CG algorithm for a V-polytope (3d simplex)
//' VP = gen_simplex(3,'V')
//' vol = volume(VP, settings = list("algorithm" = "CG"))
//'
//' # calling CG algorithm for a 2-dimensional zonotope defined as the Minkowski sum of 4 segments
//' Z = gen_rand_zonotope(2, 4)
//' vol = volume(Z, settings = list("random_walk" = "RDHR", "walk_length" = 2))
//'
//' @export
// [[Rcpp::export]]
double volume (Rcpp::Reference P,
               Rcpp::Nullable<Rcpp::List> settings = R_NilValue,
               Rcpp::Nullable<std::string> rounding = R_NilValue,
               Rcpp::Nullable<double> seed = R_NilValue) {

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef IntersectionOfVpoly<Vpolytope, RNGType> InterVP;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    unsigned int n = P.field("dimension"), walkL, type = P.field("type");

    RNGType rng(n);
    if (seed.isNotNull()) {
        unsigned seed_rcpp = Rcpp::as<double>(seed);
        rng.set_seed(seed_rcpp);
    }

    bool round = false, hpoly = false;
    unsigned int win_len = 250;

    random_walks walk;
    volume_algorithms algo;
    rounding_type rounding_method;

    if (!rounding.isNotNull()) {
        rounding_method = (type == 2) ? min_ellipsoid : none;
    } else if (Rcpp::as<std::string>(rounding).compare(std::string("min_ellipsoid")) == 0) {
        rounding_method = min_ellipsoid;
    } else if (Rcpp::as<std::string>(rounding).compare(std::string("max_ellipsoid")) == 0) {
        if (type != 1) throw Rcpp::exception("This rounding method can be used only for H-polytopes!");
        rounding_method = max_ellipsoid;
    } else if (Rcpp::as<std::string>(rounding).compare(std::string("svd")) == 0) {
        rounding_method = svd;
    } else if (Rcpp::as<std::string>(rounding).compare(std::string("none")) == 0) {
        rounding_method = none;
    }

    NT e;

    if (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("algorithm")) {
        algo = CB;
        walkL = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("walk_length")) ? 1 : Rcpp::as<int>(
                Rcpp::as<Rcpp::List>(settings)["walk_length"]);
        e = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("error")) ? 0.1 : Rcpp::as<NT>(
                Rcpp::as<Rcpp::List>(settings)["error"]);
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["algorithm"]).compare(std::string("SOB")) == 0) {
        algo = SOB;
        walkL = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("walk_length")) ? 10 + n / 10 : Rcpp::as<int>(
                Rcpp::as<Rcpp::List>(settings)["walk_length"]);
        e = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("error")) ? 1.0 : Rcpp::as<NT>(
                Rcpp::as<Rcpp::List>(settings)["error"]);
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["algorithm"]).compare(std::string("CG")) == 0) {
        algo = CG;
        walkL = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("walk_length")) ? 1 : Rcpp::as<int>(
                Rcpp::as<Rcpp::List>(settings)["walk_length"]);
        e = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("error")) ? 0.1 : Rcpp::as<NT>(
                Rcpp::as<Rcpp::List>(settings)["error"]);
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["algorithm"]).compare(std::string("CB")) == 0) {
        algo = CB;
        walkL = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("walk_length")) ? 1 : Rcpp::as<int>(
                Rcpp::as<Rcpp::List>(settings)["walk_length"]);
        e = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("error")) ? 0.1 : Rcpp::as<NT>(
                Rcpp::as<Rcpp::List>(settings)["error"]);
    } else {
        throw Rcpp::exception("Unknown method!");
    }

    if (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("random_walk")) {
        if (algo == CB) {
            walk = (type == 1) ? accelarated_billiard : billiard;
        } else {
            win_len = 4 * n * n + 500;
            if (type == 1) {
                walk = cdhr;
            } else {
                walk = rdhr;
            }
        }
    }else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["random_walk"]).compare(std::string("CDHR")) == 0) {
        walk = cdhr;
        if (algo == CB) win_len = 3*n*n+400;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["random_walk"]).compare(std::string("RDHR")) == 0) {
        walk = rdhr;
        if (algo == CB) win_len = 3*n*n+400;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["random_walk"]).compare(std::string("BaW")) == 0) {
        walk = ball_walk;
        if (algo == CB) win_len = 3*n*n+400;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["random_walk"]).compare(std::string("BiW")) == 0) {
        if (algo == CG) {
            if (type !=1){
                Rcpp::Rcout << "Billiard walk is not supported for CG algorithm. RDHR is used."<<std::endl;
                walk = rdhr;
            } else {
                Rcpp::Rcout << "Billiard walk is not supported for CG algorithm. CDHR is used."<<std::endl;
                walk = cdhr;
            }
        } else {
            walk = (type == 1) ? accelarated_billiard : billiard;
        }
    }else {
        throw Rcpp::exception("Unknown walk type!");
    }

    if (e <= 0.0) {
        throw Rcpp::exception("The error parameter has to be a positive number!");
    }

    if (Rcpp::as<Rcpp::List>(settings).containsElementNamed("win_len")) {
        win_len = Rcpp::as<int>(Rcpp::as<Rcpp::List>(settings)["win_len"]);
        if (algo == SOB) Rf_warning("input 'win_len' can be used only for CG or CB algorithms.");
    }

    switch(type) {
        case 1: {
            // Hpolytope
            Hpolytope HP(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            return generic_volume(HP, rng, walkL, e, algo, win_len, rounding_method, walk);
        }
        case 2: {
            // Vpolytope
            Vpolytope VP(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            return generic_volume(VP, rng, walkL, e, algo, win_len, rounding_method, walk);
        }
        case 3: {
            // Zonotope
            zonotope ZP(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            if (Rcpp::as<Rcpp::List>(settings).containsElementNamed("hpoly")) {
                hpoly = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(settings)["hpoly"]);
                if (hpoly && (algo == CG || algo == SOB))
                    Rf_warning("flag 'hpoly' can be used to only in MMC of CB algorithm for zonotopes.");
            } else if (ZP.num_of_generators() / ZP.dimension() < 5) {
                hpoly = true;
            } else {
                hpoly = false;
            }
            if (hpoly && algo == CB) {
                switch (walk)
                {
                case cdhr:
                    return volume_cooling_hpoly<CDHRWalk, Hpolytope>(ZP, rng, e, walkL, win_len);
                    break;
                case rdhr:
                    return volume_cooling_hpoly<RDHRWalk, Hpolytope>(ZP, rng, e, walkL, win_len);
                    break;
                case ball_walk:
                    return volume_cooling_hpoly<BallWalk, Hpolytope>(ZP, rng, e, walkL, win_len);
                    break;
                case billiard:
                    return volume_cooling_hpoly<BilliardWalk, Hpolytope>(ZP, rng, e, walkL, win_len);
                    break;
                case accelarated_billiard:
                    return volume_cooling_hpoly<AcceleratedBilliardWalk, Hpolytope>(ZP, rng, e, walkL, win_len);
                    break;
                default:
                    throw Rcpp::exception("This random walk can not be used by CB algorithm!");
                    break;
                }
            }
            return generic_volume(ZP, rng, walkL, e, algo, win_len, rounding_method, walk);
        }
        case 4: {
            // Intersection of two V-polytopes
            Vpolytope VP1(n, Rcpp::as<MT>(P.field("V1")), VT::Ones(Rcpp::as<MT>(P.field("V1")).rows()));
            Vpolytope VP2(n, Rcpp::as<MT>(P.field("V2")), VT::Ones(Rcpp::as<MT>(P.field("V2")).rows()));
            InterVP VPcVP;
            if (!seed.isNotNull()) {
                InterVP VPcVP(VP1, VP2);
            } else {
                unsigned seed3 = Rcpp::as<double>(seed);
                InterVP VPcVP(VP1, VP2, seed3);
            }
            if (!VPcVP.is_feasible()) throw Rcpp::exception("Empty set!");
            return generic_volume(VPcVP, rng, walkL, e, algo, win_len, rounding_method, walk);
        }
    }

    return 0;
}

