// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.


#ifndef SAMPLE_ONLY_H
#define SAMPLE_ONLY_H

template <typename WalkTypePolicy,
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename Point
        >
void uniform_sampling(PointList &randPoints,
                   Polytope &P,
                   RandomNumberGenerator &rng,
                   const unsigned int &walk_len,
                   const unsigned int &rnum,
                   const Point &starting_point,
                   unsigned int const& nburns)
{

    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    //RandomNumberGenerator rng(P.dimension());
    PushBackWalkPolicy push_back_policy;

    Point p = starting_point;

    typedef RandomPointGenerator <walk> RandomPointGenerator;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, nburns, walk_len, randPoints,
                                    push_back_policy, rng);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, rnum, walk_len, randPoints,
                                push_back_policy, rng);


}

template <
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename WalkTypePolicy,
        typename Point
>
void uniform_sampling(PointList &randPoints,
                   Polytope &P,
                   RandomNumberGenerator &rng,
                   WalkTypePolicy &WalkType,
                   const unsigned int &walk_len,
                   const unsigned int &rnum,
                   const Point &starting_point,
                   unsigned int const& nburns)
{
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    //RandomNumberGenerator rng(P.dimension());
    PushBackWalkPolicy push_back_policy;
    typedef RandomPointGenerator<walk> RandomPointGenerator;

    Point p = starting_point;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, nburns, walk_len, randPoints,
                                    push_back_policy, rng, WalkType.param);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, rnum, walk_len, randPoints,
                                push_back_policy, rng, WalkType.param);
}


template
<
        typename WalkTypePolicy,
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename Point
>
void uniform_sampling_boundary(PointList &randPoints,
                            Polytope &P,
                            RandomNumberGenerator &rng,
                            const unsigned int &walk_len,
                            const unsigned int &rnum,
                            const Point &starting_point,
                            unsigned int const& nburns)
{
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    //RandomNumberGenerator rng(P.dimension());
    PushBackWalkPolicy push_back_policy;

    Point p = starting_point;

    typedef BoundaryRandomPointGenerator <walk> BoundaryRandomPointGenerator;
    if (nburns > 0) {
        BoundaryRandomPointGenerator::apply(P, p, nburns, walk_len,
                                            randPoints, push_back_policy, rng);
        randPoints.clear();
    }
    unsigned int n = rnum / 2;
    BoundaryRandomPointGenerator::apply(P, p, rnum / 2, walk_len,
                                        randPoints, push_back_policy, rng);

}


template
<
        typename WalkTypePolicy,
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename NT,
        typename Point
>
void gaussian_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const NT &a,
                       const Point &starting_point,
                       unsigned int const& nburns)
{

    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    //RandomNumberGenerator rng(P.dimension());
    PushBackWalkPolicy push_back_policy;

    Point p = starting_point;

    typedef GaussianRandomPointGenerator <walk> RandomPointGenerator;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, a, nburns, walk_len, randPoints,
                                    push_back_policy, rng);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, a, rnum, walk_len, randPoints,
                                push_back_policy, rng);


}


template <
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename WalkTypePolicy,
        typename NT,
        typename Point
        >
void gaussian_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       WalkTypePolicy &WalkType,
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const NT &a,
                       const Point &starting_point,
                       unsigned int const& nburns)
{

    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    //RandomNumberGenerator rng(P.dimension());
    PushBackWalkPolicy push_back_policy;

    Point p = starting_point;

    typedef GaussianRandomPointGenerator <walk> RandomPointGenerator;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, a, nburns, walk_len, randPoints,
                                    push_back_policy, rng, WalkType.param);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, a, rnum, walk_len, randPoints,
                                push_back_policy, rng, WalkType.param);
}

template <
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename WalkTypePolicy,
        typename NT,
        typename Point,
        typename NegativeGradientFunctor,
        typename NegativeLogprobFunctor,
        typename Solver
        >
void logconcave_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const Point &starting_point,
                       unsigned int const& nburns,
                       NegativeGradientFunctor &F,
                       NegativeLogprobFunctor &f)
{
    typedef typename WalkTypePolicy::template Walk
            <
                    Point,
                    Polytope,
                    RandomNumberGenerator,
                    NegativeGradientFunctor,
                    NegativeLogprobFunctor,
                    Solver
            > walk;

    typedef typename WalkTypePolicy::template parameters
            <
                    NT,
                    NegativeGradientFunctor
            > walk_params;

    // Initialize random walk parameters
    unsigned int dim = starting_point.dimension();
    walk_params params(F, dim);

    if (F.params.eta > 0) {
        params.eta = F.params.eta;
    }

    PushBackWalkPolicy push_back_policy;

    Point p = starting_point;

    typedef LogconcaveRandomPointGenerator<walk> RandomPointGenerator;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, nburns, walk_len, randPoints,
                                    push_back_policy, rng, F, f, params);

        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, rnum, walk_len, randPoints,
                                push_back_policy, rng, F, f, params);

}


template
<
        typename WalkTypePolicy,
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename NT,
        typename Point
>
void exponential_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const Point &c,
                       const NT &a,
                       const Point &starting_point,
                       unsigned int const& nburns)
{

    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    PushBackWalkPolicy push_back_policy;

    Point p = starting_point;

    typedef ExponentialRandomPointGenerator <walk> RandomPointGenerator;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, c, a, nburns, walk_len, randPoints,
                                    push_back_policy, rng);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, c, a, rnum, walk_len, randPoints,
                                push_back_policy, rng);
}


template <
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename WalkTypePolicy,
        typename NT,
        typename Point
        >
void exponential_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       WalkTypePolicy &WalkType,
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const Point &c,
                       const NT &a,
                       const Point &starting_point,
                       unsigned int const& nburns)
{

    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    PushBackWalkPolicy push_back_policy;

    Point p = starting_point;

    typedef ExponentialRandomPointGenerator <walk> RandomPointGenerator;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, c, a, nburns, walk_len, randPoints,
                                    push_back_policy, rng, WalkType.param);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, c, a, rnum, walk_len, randPoints,
                                push_back_policy, rng, WalkType.param);
}


#endif
