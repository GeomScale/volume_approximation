// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_ACCELERATED_IMPROVED_BILLIARD_WALK_HPP
#define RANDOM_WALKS_ACCELERATED_IMPROVED_BILLIARD_WALK_HPP

#include "sampling/sphere.hpp"


// Billiard walk which accelarates each step for uniform distribution

struct AcceleratedBilliardWalk
{
    AcceleratedBilliardWalk(double L)
            :   param(L, true)
    {}

    AcceleratedBilliardWalk()
            :   param(0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set)
                :   m_L(L), set_L(set)
        {}
        double m_L;
        bool set_L;
    };

    struct update_parameters
    {
        update_parameters()
                :   facet_prev(0), hit_ball(false), inner_vi_ak(0.0), ball_inner_norm(0.0)
        {}
        int facet_prev;
        bool hit_ball;
        double inner_vi_ak;
        double ball_inner_norm;
    };

    parameters param;


    template
            <
                    typename Polytope,
                    typename RandomNumberGenerator
            >
    struct Walk
    {
        typedef typename Polytope::PointType Point;
        typedef typename Polytope::MT MT;
        typedef typename Point::FT NT;

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, Point const& p, RandomNumberGenerator &rng)
        {
            _update_parameters = update_parameters();
            _L = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            initialize(P, p, rng);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, Point const& p, RandomNumberGenerator &rng,
             parameters const& params)
        {
            _update_parameters = update_parameters();
            _L = params.set_L ? params.m_L
                              : compute_diameter<GenericPolytope>
                                ::template compute<NT>(P);
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            initialize(P, p, rng);
        }

        template
                <
                        typename GenericPolytope
                >
        inline void apply(GenericPolytope const& P,
                          Point &p,   // a point to start
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            NT T;
            const NT dl = 0.995;
            int it;

            for (auto j=0u; j<walk_length; ++j)
            {
                T = -std::log(rng.sample_urdist()) * _L;
                _v = GetDirection<Point>::apply(n, rng);
                Point p0 = _p;

                it = 0;
                std::pair<NT, int> pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _update_parameters);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    continue;
                }

                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, _update_parameters);
                it++;

                while (it < 100*n)
                {
                    std::pair<NT, int> pbpair
                            = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_parameters);
                    if (T <= pbpair.first) {
                        _p += (T * _v);
                        _lambda_prev = T;
                        break;
                    }
                    _lambda_prev = dl * pbpair.first;
                    _p += (_lambda_prev * _v);
                    T -= _lambda_prev;
                    P.compute_reflection(_v, _p, _update_parameters);
                    it++;
                }
                if (it == 100*n) _p = p0;
            }
            p = _p;
        }


        template
        <
            typename GenericPolytope
        >
        inline void get_starting_point(GenericPolytope const& P,
                           Point const& center,
                           Point &q,
                           unsigned int const& walk_length,
                           RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            NT T, Lmax = _L, max_dist, rad = 0.0, radius = P.InnerBall().second;

            q = GetPointInDsphere<Point>::apply(n, radius, rng);
            q += center;
            initialize(P, q, rng);

            apply(P, q, walk_length, rng);
        }


        template
        <
            typename GenericPolytope
        >
        inline void parameters_burnin(GenericPolytope const& P, 
                                     Point const& center,
                                     unsigned int const& num_points,
                                     unsigned int const& walk_length,
                                     RandomNumberGenerator &rng)
        {
            Point p = center;
            std::vector<Point> pointset;
            pointset.push_back(center);
            pointset.push_back(_p);
            NT rad = NT(0), max_dist, Lmax;

            for (int i = 0; i < num_points; i++) 
            {
                apply(P, p, walk_length, rng);
                max_dist = get_max_distance(pointset, p, rad);
                if (max_dist > Lmax) 
                {
                    Lmax = max_dist;
                }
                pointset.push_back(p);
            }

            if (Lmax > _L) {
                if (P.dimension() <= 2500) 
                {
                    update_delta(Lmax);
                }
            }
            pointset.clear();
        }

        

        inline void update_delta(NT L)
        {
            _L = L;
        }

    private :

        template
                <
                        typename GenericPolytope
                >
        inline void initialize(GenericPolytope const& P,
                               Point const& p,
                               RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            const NT dl = 0.995;
            _lambdas.setZero(P.num_of_hyperplanes());
            _Av.setZero(P.num_of_hyperplanes());
            _p = p;
            _v = GetDirection<Point>::apply(n, rng);

            NT T = -std::log(rng.sample_urdist()) * _L;
            Point p0 = _p;
            int it = 0;

            std::pair<NT, int> pbpair
                    = P.line_first_positive_intersect(_p, _v, _lambdas, _Av, _update_parameters);
            if (T <= pbpair.first) {
                _p += (T * _v);
                _lambda_prev = T;
                return;
            }
            _lambda_prev = dl * pbpair.first;
            _p += (_lambda_prev * _v);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, _update_parameters);

            while (it <= 100*n)
            {
                std::pair<NT, int> pbpair
                        = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_parameters);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    break;
                } else if (it == 100*n) {
                    _lambda_prev = rng.sample_urdist() * pbpair.first;
                    _p += (_lambda_prev * _v);
                    break;
                }
                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, _update_parameters);
                it++;
            }
        }

        inline double get_max_distance(std::vector<Point> &pointset, Point const& q, double &rad) 
        {
            double dis = -1.0, temp_dis;
            int jj = 0;
            for (auto vecit = pointset.begin(); vecit!=pointset.end(); vecit++, jj++) 
            {
                temp_dis = (q.getCoefficients() - (*vecit).getCoefficients()).norm();
                if (temp_dis > dis) {
                    dis = temp_dis;
                }
                if (jj == 0) {
                    if (temp_dis > rad) {
                        rad = temp_dis;
                    }
                }
            }
            return dis;
        }

        double _L;
        Point _p;
        Point _v;
        NT _lambda_prev;
        MT _AA;
        update_parameters _update_parameters;
        typename Point::Coeff _lambdas;
        typename Point::Coeff _Av;
    };

};


#endif


