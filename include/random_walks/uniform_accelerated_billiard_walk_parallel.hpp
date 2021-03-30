// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_ACCELERATED_IMPROVED_BILLIARD_WALK_PARALLEL_HPP
#define RANDOM_WALKS_ACCELERATED_IMPROVED_BILLIARD_WALK_PARALLEL_HPP

#include "sampling/sphere.hpp"


// Billiard walk which accelarates each step for uniform distribution and can be used for a parallel use by threads

struct AcceleratedBilliardWalkParallel
{
    AcceleratedBilliardWalkParallel()
    {}


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

    template<typename NT, typename Point>
    struct thread_parameters
    {
        thread_parameters(unsigned int d, unsigned int m)
        {
            update_step_parameters = update_parameters();
            p = Point(d);
            v = Point(d);
            lambdas.setZero(m);
            Av.setZero(m);
            lambda_prev = NT(0);
        }

        update_parameters update_step_parameters;
        Point p;
        Point v;
        NT lambda_prev;
        typename Point::Coeff lambdas;
        typename Point::Coeff Av;
    };


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
        Walk(GenericPolytope const& P)
        {
            _L = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            _p0 = Point(P.dimension());
            _rho = 1000 * P.dimension();
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, NT const& L)
        {
            _L = L > NT(0) ? L
                              : compute_diameter<GenericPolytope>
                                ::template compute<NT>(P);
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            _p0 = Point(P.dimension());
            _rho = 1000 * P.dimension();
        }

        template
        <
            typename GenericPolytope,
            typename thread_params
        >
        inline void apply(GenericPolytope const& P,
                          thread_params &params,   // a point to start
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
                params.v = GetDirection<Point>::apply(n, rng);
                _p0 = params.p;

                it = 0;
                std::pair<NT, int> pbpair = P.line_positive_intersect(params.p, params.v, params.lambdas, params.Av, 
                                                                      params.lambda_prev, params.update_step_parameters);
                if (T <= pbpair.first) 
                {
                    params.p += (T * params.v);
                    params.lambda_prev = T;
                    continue;
                }

                params.lambda_prev = dl * pbpair.first;
                params.p += (params.lambda_prev * params.v);
                T -= params.lambda_prev;
                P.compute_reflection(params.v, params.p, params.update_step_parameters);
                it++;

                while (it < _rho)
                {
                    std::pair<NT, int> pbpair
                            = P.line_positive_intersect(params.p, params.v, params.lambdas, params.Av, 
                                                        params.lambda_prev, _AA, params.update_step_parameters);
                    if (T <= pbpair.first) {
                        params.p += (T * params.v);
                        params.lambda_prev = T;
                        break;
                    }
                    params.lambda_prev = dl * pbpair.first;
                    params.p += (params.lambda_prev * params.v);
                    T -= params.lambda_prev;
                    P.compute_reflection(params.v, params.p, params.update_step_parameters);
                    it++;
                }
                if (it == _rho) params.p = _p0;
            }
        }


        template
        <
            typename GenericPolytope,
            typename thread_params
        >
        inline void get_starting_point(GenericPolytope const& P,
                                       Point const& center,
                                       thread_params &params,
                                       unsigned int const& walk_length,
                                       RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            NT T, Lmax = _L, max_dist, rad = 0.0, radius = P.InnerBall().second;

            params.p = GetPointInDsphere<Point>::apply(n, radius, rng);
            params.p += center;
            initialize(P, params, rng);

            apply(P, params, walk_length, rng);
        }


        template
        <
            typename GenericPolytope,
            typename thread_params
        >
        inline void parameters_burnin(GenericPolytope const& P, 
                                      Point const& center,
                                      unsigned int const& num_points,
                                      unsigned int const& walk_length,
                                      RandomNumberGenerator &rng,
                                      thread_params &params)
        {
            std::vector<Point> pointset;
            pointset.push_back(center);
            params.p = center;
            initialize(P, params, rng);
            pointset.push_back(params.p);
            NT rad = NT(0), max_dist, Lmax;

            for (int i = 0; i < num_points; i++) 
            {
                apply(P, params, walk_length, rng);
                max_dist = get_max_distance(pointset, params.p, rad);
                if (max_dist > Lmax) 
                {
                    Lmax = max_dist;
                }
                pointset.push_back(params.p);
            }

            if (Lmax > _L) 
            {
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
            typename GenericPolytope,
            typename thread_params
        >
        inline void initialize(GenericPolytope const& P,
                               thread_params &params,
                               RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            const NT dl = 0.995;
            params.v = GetDirection<Point>::apply(n, rng);

            NT T = -std::log(rng.sample_urdist()) * _L;
            int it = 0;

            std::pair<NT, int> pbpair
                    = P.line_first_positive_intersect(params.p, params.v, params.lambdas, 
                                                      params.Av, params.update_step_parameters);
            if (T <= pbpair.first) {
                params.p += (T * params.v);
                params.lambda_prev = T;
                return;
            }
            params.lambda_prev = dl * pbpair.first;
            params.p += (params.lambda_prev * params.v);
            T -= params.lambda_prev;
            P.compute_reflection(params.v, params.p, params.update_step_parameters);

            while (it <= _rho)
            {
                std::pair<NT, int> pbpair
                        = P.line_positive_intersect(params.p, params.v, params.lambdas, params.Av, 
                                                    params.lambda_prev, _AA, params.update_step_parameters);
                if (T <= pbpair.first) {
                    params.p += (T * params.v);
                    params.lambda_prev = T;
                    break;
                } else if (it == _rho) {
                    params.lambda_prev = rng.sample_urdist() * pbpair.first;
                    params.p += (params.lambda_prev * params.v);
                    break;
                }
                params.lambda_prev = dl * pbpair.first;
                params.p += (params.lambda_prev * params.v);
                T -= params.lambda_prev;
                P.compute_reflection(params.v, params.p, params.update_step_parameters);
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

        NT _L;
        MT _AA;
        Point _p0;
        unsigned int _rho;
    };

};


#endif


