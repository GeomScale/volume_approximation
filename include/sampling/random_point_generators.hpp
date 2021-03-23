// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SAMPLERS_RANDOM_POINT_GENERATORS_HPP
#define SAMPLERS_RANDOM_POINT_GENERATORS_HPP

template
<
    typename Walk
>
struct RandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator,
        typename Parameters
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        Walk walk(P, p, rng, parameters);
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }

    template
    <
            typename Polytope,
            typename Point,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, rng);
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }
};

template
<
    typename Walk
>
struct GaussianRandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename NT,
        typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator
    >
    static void apply(Polytope const& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, a_i, rng);
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, a_i, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }

    template
    <
            typename Polytope,
            typename Point,
            typename NT,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator,
            typename Parameters
    >
    static void apply(Polytope const& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        Walk walk(P, p, a_i, rng, parameters);

        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, a_i, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }
};



template <typename Walk>
struct BoundaryRandomPointGenerator
{
    template
    <
            typename Polytope,
            typename Point,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator
    >
    static void apply(Polytope const& P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, rng);
        Point p1(P.dimension()), p2(P.dimension());
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p1, p2, walk_length, rng);
            policy.apply(randPoints, p1);
            policy.apply(randPoints, p2);
        }
    }
};


template
<
    typename Walk
>
struct LogconcaveRandomPointGenerator
{

    template
    <
            typename Polytope,
            typename Point,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator,
            typename NegativeGradientFunctor,
            typename NegativeLogprobFunctor,
            typename Parameters
    >
    static void apply(Polytope &P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      NegativeGradientFunctor &F,
                      NegativeLogprobFunctor &f,
                      Parameters &parameters)
    {
        typedef double NT;
        Walk walk(&P, p, F, f, parameters);

        for (unsigned int i = 0; i < rnum; ++i)
        {
            // Gather one sample
            walk.apply(rng, walk_length);

            // Use PushBackWalkPolicy
            policy.apply(randPoints, walk.x);
        }
    }
};


template
<
    typename Walk
>
struct ExponentialRandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename NT,
        typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator
    >
    static void apply(Polytope const& P,
                      Point &p,   // a point to start
                      Point const& c,   // bias function
                      NT const& T, // temperature/variance
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, c, T, rng);
        bool success;
        for (unsigned int i=0; i<rnum; ++i)
        {
            success = walk.template apply(P, p, walk_length, rng);
            if (!success) {
                return;
            }
            policy.apply(randPoints, p);
        }
    }

    template
    <
            typename Polytope,
            typename Point,
            typename NT,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator,
            typename Parameters
    >
    static void apply(Polytope const& P,
                      Point &p,   // a point to start
                      Point const& c,   // bias function
                      NT const& T, // temperature/variance
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        Walk walk(P, p, c, T, rng, parameters);
        bool success;

        for (unsigned int i=0; i<rnum; ++i)
        {
            success = walk.template apply(P, p, walk_length, rng);
            if (!success) {
                return;
            }
            policy.apply(randPoints, p);
        }
    }
};



#endif // SAMPLERS_RANDOM_POINT_GENERATORS_HPP
