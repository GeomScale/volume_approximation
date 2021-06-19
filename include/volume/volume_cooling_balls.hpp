// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLUME_COOLING_BALLS_HPP
#define VOLUME_COOLING_BALLS_HPP

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/erf.hpp>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#ifndef VOLESTIPY
    #include "convex_bodies/vpolytope.h"
    #include "convex_bodies/zpolytope.h"
    #include "convex_bodies/vpolyintersectvpoly.h"
#endif
#include "convex_bodies/ballintersectconvex.h"
#include "sampling/random_point_generators.hpp"


////////////////////////////////////
// ball annealing

template <typename NT>
struct cooling_ball_parameters
{
    cooling_ball_parameters(unsigned int win_len)
        :   lb(0.1)
        ,   ub(0.15)
        ,   p(0.75)
        ,   rmax(0)
        ,   alpha(0.2)
        ,   win_len(win_len)
        ,   N(125)
        ,   nu(10)
        ,   window2(false)
    {}

    NT lb;
    NT ub;
    NT p;
    NT rmax;
    NT alpha;
    int win_len;
    int N;
    int nu;
    bool window2;
};

/// Helpers

template <typename Point, typename ConvexBody, typename PointList, typename NT>
bool check_convergence(ConvexBody const& P,
                       PointList const& randPoints,
                       bool& too_few,
                       NT& ratio,
                       int const& nu,
                       bool const& precheck,
                       bool const& lastball,
                       cooling_ball_parameters<NT> const& parameters)
{
    NT alpha = parameters.alpha;
    std::vector<NT> ratios;
    std::pair<NT,NT> mv;
    int m = randPoints.size()/nu;
    int i = 1;
    NT T;
    NT rs;
    NT alpha_check = 0.01;
    size_t countsIn = 0;

    for (auto pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++)
    {
        if (P.is_in(*pit)==-1) countsIn++;
        if (i % m == 0)
        {
            ratios.push_back(NT(countsIn)/m);
            countsIn = 0;
            if (ratios.size()>1 && precheck)
            {
                boost::math::students_t dist(ratios.size() - 1);
                mv = get_mean_variance(ratios);
                ratio = mv.first;
                rs = std::sqrt(mv.second);
                T = rs * (boost::math::quantile
                            (boost::math::complement(dist, alpha_check / 2.0))
                          / std::sqrt(NT(ratios.size())));
                if (ratio + T < parameters.lb)
                {
                    too_few = true;
                    return false;
                } else if (ratio - T > parameters.ub) return false;
            }
        }
    }

    if (precheck) alpha *= 0.5;
    mv = get_mean_variance(ratios);
    ratio = mv.first;
    rs = std::sqrt(mv.second);
    boost::math::students_t dist(nu - 1);
    T = rs * (boost::math::quantile(boost::math::complement(dist, alpha))
           / std::sqrt(NT(nu)));
    if (ratio > parameters.lb + T)
    {
        if (lastball) return true;
        if ((precheck && ratio < parameters.ub - T)
        || (!precheck && ratio < parameters.ub + T)) return true;
        return false;
    }
    too_few = true;
    return false;
}


template <typename Polytope, typename Ball, typename NT, typename RNG>
bool get_first_ball(Polytope const& P,
                    Ball& B0,
                    NT& ratio,
                    NT const& radius_input,
                    cooling_ball_parameters<NT> const& parameters,
                    RNG& rng) {
    const unsigned max_iterarions = 20;
    NT tolerance = 0.00000000001;
    typedef typename Polytope::PointType Point;
    int n = P.dimension();
    int iter = 1;
    bool bisection_int = false;
    bool pass = false;
    bool too_few = false;
    std::list <Point> randPoints;
    NT rmax = parameters.rmax;
    NT sqrt_n = std::sqrt(NT(n));
    NT radius1 = radius_input;

    if (rmax > 0.0) {
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(GetPointInDsphere<Point>::apply(n, rmax, rng));
        }
        pass = check_convergence<Point>(P, randPoints, too_few, ratio,
                                        10, true, false, parameters);
        if (pass || !too_few) {
            B0 = Ball(Point(n), rmax * rmax);
            return true;
        }
        bisection_int = true;
    } else {
        rmax = 2 * sqrt_n * radius1;
    }
    NT radius = radius1;

    while (!bisection_int) {
        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(GetPointInDsphere<Point>::apply(n, rmax, rng));
        }

        if (check_convergence<Point>(P, randPoints, too_few, ratio, 10,
                                     true, false, parameters)) {
            B0 = Ball(Point(n), rmax * rmax);
            return true;
        }

        if (too_few) break;
        radius1 = rmax;
        rmax = rmax + 2 * sqrt_n * radius;
    }

    NT rad_med;
    NT rad0 = radius1;
    NT rad_m = rmax;

    while (iter <= max_iterarions) {
        rad_med = 0.5 * (radius1 + rmax);
        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(GetPointInDsphere<Point>::apply(n, rad_med, rng));
        }

        if (check_convergence<Point>(P, randPoints, too_few, ratio, 10,
                                     true, false, parameters)) {
            B0 = Ball(Point(n), rad_med * rad_med);
            return true;
        }

        if (too_few) {
            rmax = rad_med;
        } else {
            radius1 = rad_med;
        }

        if (rmax - radius1 < tolerance) {
            radius1 = rad0;
            rmax = rad_m;
            iter++;
        }
    }
    return false;
}


template <typename Point, typename ball, typename PointList, typename NT>
bool get_next_zonotopeball(std::vector<ball>& BallSet,
                           PointList const& randPoints,
                           NT const& rad_min,
                           std::vector<NT>& ratios,
                           cooling_ball_parameters<NT> const& parameters)
{
    const unsigned max_iterarions = 20;
    NT tolerance = 0.00000000001;
    int n = (*randPoints.begin()).dimension();
    int iter = 1;
    bool too_few;
    NT radmax = NT(0);
    NT radmin = rad_min;

    for (auto rpit = randPoints.begin();  rpit!=randPoints.end(); ++rpit)
    {
        NT pnorm = (*rpit).squared_length();
        if (pnorm > radmax) radmax = pnorm;
    }
    ball Biter;
    radmax = std::sqrt(radmax);
    NT radmin_init = radmin;
    NT radmax_init = radmax;

    while (iter <= max_iterarions)
    {
        NT rad = 0.5 * (radmin + radmax);
        Biter = ball(Point(n), rad * rad);
        too_few = false;

        NT ratio;
        if (check_convergence<Point>(Biter, randPoints, too_few, ratio,
                                     parameters.nu, false, false, parameters))
        {
            BallSet.push_back(Biter);
            ratios.push_back(ratio);
            return true;
        }

        if (too_few)
        {
            radmin = rad;
        } else
        {
            radmax = rad;
        }

        if (radmax-radmin < tolerance)
        {
            radmin = radmin_init;
            radmax = radmax_init;
            iter++;
        }
    }
    return false;
}


template
<
    typename RandomPointGenerator,
    typename PolyBall,
    typename ball,
    typename Polytope,
    typename NT,
    typename RNG
>
bool get_sequence_of_polytopeballs(Polytope& P,
                                   std::vector<ball>& BallSet,
                                   std::vector<NT>& ratios,
                                   int const& Ntot,
                                   NT const& radius,
                                   unsigned int const& walk_length,
                                   cooling_ball_parameters<NT> const& parameters,
                                   RNG& rng)
{

    typedef typename Polytope::PointType Point;
    bool fail;
    int n = P.dimension();
    NT ratio;
    NT ratio0;
    std::list<Point> randPoints;
    ball B0;
    Point q(n);

    if ( !get_first_ball(P, B0, ratio, radius, parameters, rng) )
    {
        return false;
    }

    ratio0 = ratio;

    PushBackWalkPolicy push_back_policy;
    RandomPointGenerator::apply(P, q, Ntot, walk_length,
                                randPoints, push_back_policy, rng);

    if (check_convergence<Point>(B0, randPoints,
                                 fail, ratio, parameters.nu,
                                 false, true, parameters))
    {
        ratios.push_back(ratio);
        BallSet.push_back(B0);
        ratios.push_back(ratio0);
        return true;
    }
    if ( !get_next_zonotopeball<Point>(BallSet, randPoints, B0.radius(), ratios,
                                       parameters) )
    {
        return false;
    }

    while (true)
    {
        PolyBall zb_it(P, BallSet[BallSet.size()-1]);
        q.set_to_origin();
        randPoints.clear();

        RandomPointGenerator::apply(zb_it, q, Ntot, walk_length,
                                    randPoints, push_back_policy, rng);
        if (check_convergence<Point>(B0, randPoints, fail, ratio, parameters.nu,
                                     false, true, parameters))
        {
            ratios.push_back(ratio);
            BallSet.push_back(B0);
            ratios.push_back(ratio0);
            return true;
        }
        if ( !get_next_zonotopeball<Point>(BallSet, randPoints, B0.radius(),
                                           ratios, parameters) )
        {
            return false;
        }
    }
}


////////////////////////////////////
///
/// ratio estimation

template <typename NT>
bool is_max_error(NT const& a, NT const& b, NT const& error)
{
    return ((b-a)/a<error/2.0) ? true : false;
}

template <typename NT>
struct estimate_ratio_parameters
{
public:

    estimate_ratio_parameters(unsigned int W_len, unsigned int N, NT ratio)
        :   min_val(std::numeric_limits<NT>::lowest())
        ,   max_val(std::numeric_limits<NT>::max())
        ,   max_iterations_estimation(10000000)
        ,   min_index(W_len-1)
        ,   max_index(W_len-1)
        ,   W(W_len)
        ,   index(0)
        ,   tot_count(N)
        ,   count_in(N * ratio)
        ,   iter(0)
        ,   last_W(std::vector<NT>(W_len))
        ,   minmaxIt(last_W.begin())
    {}

    NT min_val;
    NT max_val;
    const unsigned int max_iterations_estimation;
    unsigned int min_index;
    unsigned int max_index;
    unsigned int W;
    unsigned int index;
    size_t tot_count;
    size_t count_in;
    unsigned int iter;
    std::vector<NT> last_W;
    typename std::vector<NT>::iterator minmaxIt;
};

template <typename Pollyball, typename Point, typename NT>
bool estimate_ratio_generic(Pollyball const& Pb2, Point const& p, NT const& error,
       estimate_ratio_parameters<NT> &ratio_parameters)
{
    if (ratio_parameters.iter++ <= ratio_parameters.max_iterations_estimation)
    {
        if (Pb2.is_in(p) == -1) ratio_parameters.count_in = ratio_parameters.count_in + 1.0;

        ratio_parameters.tot_count = ratio_parameters.tot_count + 1.0;
        NT val = NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
        ratio_parameters.last_W[ratio_parameters.index] = val;

        if (val <= ratio_parameters.min_val)
        {
            ratio_parameters.min_val = val;
            ratio_parameters.min_index = ratio_parameters.index;
        } else if (ratio_parameters.min_index == ratio_parameters.index)
        {
            ratio_parameters.minmaxIt = std::min_element(ratio_parameters.last_W.begin(), ratio_parameters.last_W.end());
            ratio_parameters.min_val = (*ratio_parameters.minmaxIt);
            ratio_parameters.min_index = std::distance(ratio_parameters.last_W.begin(), ratio_parameters.minmaxIt);
        }

        if (val >= ratio_parameters.max_val)
        {
            ratio_parameters.max_val = val;
            ratio_parameters.max_index = ratio_parameters.index;
        } else if (ratio_parameters.max_index == ratio_parameters.index)
        {
            ratio_parameters.minmaxIt = std::max_element(ratio_parameters.last_W.begin(), ratio_parameters.last_W.end());
            ratio_parameters.max_val = (*ratio_parameters.minmaxIt);
            ratio_parameters.max_index = std::distance(ratio_parameters.last_W.begin(), ratio_parameters.minmaxIt);
        }

        if ( (ratio_parameters.max_val - ratio_parameters.min_val) / ratio_parameters.max_val <= error/2.0 )
        {
            return true;
        }

        ratio_parameters.index = ratio_parameters.index % ratio_parameters.W + 1;
        if (ratio_parameters.index == ratio_parameters.W) ratio_parameters.index = 0;

        return false;
    }
    return true;
}


template
<
        typename WalkType,
        typename Point,
        typename PolyBall1,
        typename PolyBall2,
        typename NT,
        typename RNG
>
NT estimate_ratio(PolyBall1 const& Pb1,
                  PolyBall2 const& Pb2,
                  NT const& ratio,
                  NT const& error,
                  unsigned int const& W,
                  unsigned int const& Ntot,
                  unsigned int const& walk_length,
                  RNG& rng)
{
    estimate_ratio_parameters<NT> ratio_parameters(W, Ntot, ratio);
    unsigned int n = Pb1.dimension();
    Point p(n);
    WalkType walk(Pb1, p, rng);

    do
    {
        walk.template apply(Pb1, p, walk_length, rng);
    } while(!estimate_ratio_generic(Pb2, p, error, ratio_parameters));

    return NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
}

template
<       typename Point,
        typename ball,
        typename PolyBall,
        typename NT,
        typename RNG
>
NT estimate_ratio(ball const& B,
                  PolyBall const& Pb2,
                  NT const& ratio,
                  NT const& error,
                  int const& W,
                  int const& Ntot,
                  RNG& rng)
{
    estimate_ratio_parameters<NT> ratio_parameters(W, Ntot, ratio);
    unsigned int n = B.dimension();
    Point p(n);
    NT radius = B.radius();

    do
    {
        p = GetPointInDsphere<Point>::apply(n, radius, rng);
    } while(!estimate_ratio_generic(Pb2, p, error, ratio_parameters));

    return NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
}

//--------------------------------------------------------------------------

template <typename NT>
struct estimate_ratio_interval_parameters
{
public:
    estimate_ratio_interval_parameters(unsigned int W_len,
                                       unsigned int N,
                                       NT ratio)
        :    mean(0)
        ,    sum_sq(0)
        ,    sum(0)
        ,    s(0)
        ,    max_iterations_estimation(10000000)
        ,    W(W_len)
        ,    index(0)
        ,    tot_count(N)
        ,    count_in(N * ratio)
        ,    iter(0)
        ,    last_W(std::vector<NT>(W_len))
    {}

    NT mean;
    NT sum_sq;
    NT sum;
    NT s;
    const unsigned int max_iterations_estimation;
    unsigned int W;
    unsigned int index;
    size_t tot_count;
    size_t count_in;
    unsigned int iter;
    std::vector<NT> last_W;
};

template <typename Pollyball, typename Point, typename NT>
void full_sliding_window(Pollyball const& Pb2,
                         Point const& p,
                         estimate_ratio_interval_parameters<NT>& ratio_parameters)
{
    if (Pb2.is_in(p) == -1) ratio_parameters.count_in = ratio_parameters.count_in + 1.0;

    ratio_parameters.tot_count = ratio_parameters.tot_count + 1.0;
    NT val = NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
    ratio_parameters.sum += val;
    ratio_parameters.sum_sq += val * val;
    ratio_parameters.last_W[ratio_parameters.index] = val;
    ratio_parameters.index = ratio_parameters.index % ratio_parameters.W + 1;
    if (ratio_parameters.index == ratio_parameters.W) ratio_parameters.index = 0;
}

template <typename Pollyball, typename Point, typename NT>
bool estimate_ratio_interval_generic(Pollyball const& Pb2,
                                     Point const& p,
                                     NT const& error,
                                     NT const& zp,
                                     estimate_ratio_interval_parameters
                                     <NT>& ratio_parameters)
{
    if (ratio_parameters.iter++ <= ratio_parameters.max_iterations_estimation)
    {
        if (Pb2.is_in(p) == -1) ratio_parameters.count_in = ratio_parameters.count_in + 1.0;

        ratio_parameters.tot_count = ratio_parameters.tot_count + 1.0;
        NT val = NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);

        ratio_parameters.mean = (ratio_parameters.mean
                             - ratio_parameters.last_W[ratio_parameters.index] /
                NT(ratio_parameters.W)) + val / NT(ratio_parameters.W);

        ratio_parameters.sum_sq = (ratio_parameters.sum_sq -
                ratio_parameters.last_W[ratio_parameters.index]
                * ratio_parameters.last_W[ratio_parameters.index])
                + val * val;

        ratio_parameters.sum = (ratio_parameters.sum
                                - ratio_parameters.last_W[ratio_parameters.index])
                               + val;

        ratio_parameters.s = std::sqrt((ratio_parameters.sum_sq + NT(ratio_parameters.W) *
                ratio_parameters.mean * ratio_parameters.mean - NT(2)
                                      * ratio_parameters.mean
                                      * ratio_parameters.sum) /
                                       NT(ratio_parameters.W));

        ratio_parameters.last_W[ratio_parameters.index] = val;

        ratio_parameters.index = ratio_parameters.index % ratio_parameters.W + 1;
        if (ratio_parameters.index == ratio_parameters.W)
        {
            ratio_parameters.index = 0;
        }

        if (is_max_error(val - zp * ratio_parameters.s,
                         val + zp * ratio_parameters.s,
                         error))
        {
            return true;
        }
        return false;
    }
    return true;
}

template
<
    typename Point,
    typename ball,
    typename PolyBall2,
    typename NT,
    typename RNG
>
NT estimate_ratio_interval(ball const& B,
                           PolyBall2 const& Pb2,
                           NT const& ratio,
                           NT const& error,
                           int const& W,
                           int const& Ntot,
                           NT const& prob,
                           RNG& rng)
{
    estimate_ratio_interval_parameters<NT> ratio_parameters(W, Ntot, ratio);
    boost::math::normal dist(0.0, 1.0);
    NT zp = boost::math::quantile(boost::math::complement(dist, (1.0 - prob)/2.0));
    NT radius = B.radius();

    unsigned int n = Pb2.dimension();
    Point p(n);

    for (int i = 0; i < ratio_parameters.W; ++i)
    {
        p = GetPointInDsphere<Point>::apply(n, radius, rng);
        full_sliding_window(Pb2, p, ratio_parameters);
    }

    ratio_parameters.mean = ratio_parameters.sum / NT(ratio_parameters.W);

    do {
        p = GetPointInDsphere<Point>::apply(n, radius, rng);
    } while (!estimate_ratio_interval_generic(Pb2, p, error, zp, ratio_parameters));

    return NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
}


template
<
        typename WalkType,
        typename Point,
        typename PolyBall1,
        typename PolyBall2,
        typename NT,
        typename RNG
>
NT estimate_ratio_interval(PolyBall1 const& Pb1,
                           PolyBall2 const& Pb2,
                           NT const& ratio,
                           NT const& error,
                           int const& W,
                           int const& Ntot,
                           NT const& prob,
                           unsigned int const& walk_length,
                           RNG& rng)
{
    estimate_ratio_interval_parameters<NT> ratio_parameters(W, Ntot, ratio);
    boost::math::normal dist(0.0, 1.0);
    NT zp = boost::math::quantile(boost::math::complement(dist, (1.0 - prob)/2.0));

    unsigned int n = Pb1.dimension();
    Point p(n);
    WalkType walk(Pb1, p, rng);

    for (int i = 0; i < ratio_parameters.W; ++i)
    {
        walk.template apply(Pb1, p, walk_length, rng);
        full_sliding_window(Pb2, p, ratio_parameters);
    }

    ratio_parameters.mean = ratio_parameters.sum / NT(ratio_parameters.W);

    do {
        walk.template apply(Pb1, p, walk_length, rng);
    } while (!estimate_ratio_interval_generic(Pb2, p, error, zp, ratio_parameters));

    return NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
}


template <typename NT>
NT log_gamma_function(NT x)
{
    if (x <= NT(100)) return std::log(tgamma(x));
    return (std::log(x - NT(1)) + log_gamma_function(x - NT(1)));
}


template
<
    typename WalkTypePolicy,
    typename Polytope,
    typename RandomNumberGenerator
>
std::pair<double, double> volume_cooling_balls(Polytope const& Pin,
                                       RandomNumberGenerator &rng,
                                       double const& error = 0.1,
                                       unsigned int const& walk_length = 1,
                                       unsigned int const& win_len = 300)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope <Polytope, BallType> PolyBall;
    typedef typename Polytope::VT VT;
    typedef std::list <Point> PointList;

    typedef typename WalkTypePolicy::template Walk
                                              <
                                                Polytope,
                                                RandomNumberGenerator
                                              > WalkType;
    typedef RandomPointGenerator<WalkType> RandomPointGenerator;

    auto P(Pin);
    cooling_ball_parameters<NT> parameters(win_len);

    int n = P.dimension();
    NT prob = parameters.p;
    int N_times_nu = parameters.N * parameters.nu;

    auto InnerBall = P.ComputeInnerBall();
    if (InnerBall.second < 0.0) return std::pair<NT, NT> (-1.0, 0.0);

    NT radius = InnerBall.second;
    Point c = InnerBall.first;

    std::vector<BallType> BallSet;
    std::vector<NT> ratios;

    // move the chebychev center to the origin
    // and apply the same shifting to the polytope
    P.shift(c.getCoefficients());

    if ( !get_sequence_of_polytopeballs
          <
            RandomPointGenerator,
            PolyBall
          >(P, BallSet, ratios,
            N_times_nu, radius, walk_length,
            parameters, rng) )
    {
        return std::pair<NT, NT> (-1.0, 0.0);
    }

    NT vol = (NT(n)/NT(2) * std::log(M_PI)) + NT(n)*std::log((*(BallSet.end() - 1)).radius()) - log_gamma_function(NT(n) / NT(2) + 1);

    int mm = BallSet.size() + 1;
    prob = std::pow(prob, 1.0 / NT(mm));
    NT er0 = error / (2.0 * std::sqrt(NT(mm)));
    NT er1 = (error * std::sqrt(4.0 * NT(mm) - 1)) / (2.0 * std::sqrt(NT(mm)));

    vol += (parameters.window2) ?
                std::log(estimate_ratio<Point>(*(BallSet.end() - 1),
                                               P, *(ratios.end() - 1),
                                               er0, parameters.win_len, 1200, rng))
              : std::log(estimate_ratio_interval<Point>(*(BallSet.end() - 1),
                                                        P, *(ratios.end() - 1),
                                                        er0, parameters.win_len, 1200,
                                                        prob, rng));

    auto balliter = BallSet.begin();
    auto ratioiter = ratios.begin();

    er1 = er1 / std::sqrt(NT(mm) - 1.0);

    if (*ratioiter != 1)
    {
        vol += (!parameters.window2) ?
               std::log(NT(1) / estimate_ratio_interval
                    <WalkType, Point>(P,
                                      *balliter,
                                      *ratioiter,
                                      er1,
                                      parameters.win_len,
                                      N_times_nu,
                                      prob,
                                      walk_length,
                                      rng))
            : std::log(NT(1) / estimate_ratio
                    <WalkType, Point>(P,
                                      *balliter,
                                      *ratioiter,
                                      er1,
                                      parameters.win_len,
                                      N_times_nu,
                                      walk_length,
                                      rng));
    }

    for ( ; balliter < BallSet.end() - 1; ++balliter, ++ratioiter)
    {
        PolyBall Pb(P, *balliter);
        vol += (!parameters.window2) ?
                    std::log(NT(1) / estimate_ratio_interval
                                <WalkType, Point>(Pb,
                                                  *(balliter + 1),
                                                  *(ratioiter + 1),
                                                  er1, parameters.win_len,
                                                  N_times_nu,
                                                  prob, walk_length,
                                                  rng))
                  : std::log(NT(1) / estimate_ratio
                                <WalkType, Point>(Pb,
                                                  *balliter,
                                                  *ratioiter,
                                                  er1,
                                                  parameters.win_len,
                                                  N_times_nu,
                                                  walk_length,
                                                  rng));
    }

    return std::pair<NT, NT> (vol, std::exp(vol));
}



template
<
    typename WalkTypePolicy = CDHRWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt11213b,
                                                                double>,
    typename Polytope
>
std::pair<double, double> volume_cooling_balls(Polytope const& Pin,
                                               double const& error = 0.1,
                                               unsigned int const& walk_length = 1)
{
    RandomNumberGenerator rng(Pin.dimension());
    return volume_cooling_balls<WalkTypePolicy>(Pin, rng, error, walk_length);
}


#endif // VOLUME_COOLING_BALLS_HPP
