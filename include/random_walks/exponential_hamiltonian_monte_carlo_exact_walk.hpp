// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis


// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_EXPONENTIAL_EXACT_HMC_WALK_HPP
#define RANDOM_WALKS_EXPONENTIAL_EXACT_HMC_WALK_HPP

#define IN_INSIDE_BODY_TOLLERANCE 1e-10

#include "sampling/sphere.hpp"


// Exact HMC for sampling from the Exponential distribution restricted to a convex polytope
struct ExponentialHamiltonianMonteCarloExactWalk
{
    ExponentialHamiltonianMonteCarloExactWalk(double L)
            :   param(L, true, 0, false)
    {}

    ExponentialHamiltonianMonteCarloExactWalk(double L, unsigned int rho)
            :   param(L, true, rho, true)
    {}

    ExponentialHamiltonianMonteCarloExactWalk()
            :   param(0, false, 0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set, unsigned int _rho, bool _set_rho)
                :   m_L(L), set_L(set), rho(_rho), set_rho(_set_rho)
        {}
        double m_L;
        bool set_L;
        unsigned int rho;
        bool set_rho;
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
    typedef typename Point::FT NT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p, Point const& c, NT const& T, RandomNumberGenerator &rng)
    {
        _Len = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
        _Ac = P.get_mat() * c.getCoefficients();
        _Temp = T;
        _c = c;
        _rho = 100 * P.dimension(); // upper bound for the number of reflections (experimental)
        initialize(P, p, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p, Point const& c, NT const& T, RandomNumberGenerator &rng,
         parameters const& params)
    {
        _Len = params.set_L ? params.m_L
                          : compute_diameter<GenericPolytope>
                            ::template compute<NT>(P);
        _Ac = P.get_mat() * c.getCoefficients();
        _Temp = T;
        _c = c;
        _rho = 100 * P.dimension(); // upper bound for the number of reflections (experimental)
        initialize(P, p, rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline bool apply(GenericPolytope const& P,
                      Point& p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT T;
        int failures = 0, it;
        Point p0 = _p;

        for (auto j=0u; j<walk_length; ++j)
        {
            do {
                _p = p0;
                failures++;
                if (failures == 1000) { // if the number of failures exceeds 1000 then stop
                    return false;
                }
                T = rng.sample_urdist() * _Len;
                _v = GetDirection<Point>::apply(n, rng, false);
                
                it = 0;
                while (it < _rho)
                {
                    auto pbpair = P.quadratic_positive_intersect(_p, _v, _Ac, _Temp, _lambdas,
                                                             _Av, _lambda_prev, _facet_prev);
                    if (T <= pbpair.first) {
                        _p += ((T * T) / (-2.0*_Temp)) *_c + (T * _v);
                        _lambda_prev = T;
                        break;
                    }
                    _lambda_prev = pbpair.first;
                    _p += ((_lambda_prev * _lambda_prev) / (-2.0*_Temp)) *_c + (_lambda_prev * _v);
                    T -= _lambda_prev;
                    _v += (-_lambda_prev/_Temp) * _c;
                    P.compute_reflection(_v, _p, pbpair.second);
                    it++;
                }
                
            } while (P.is_in(_p, IN_INSIDE_BODY_TOLLERANCE) == 0);
            if (it == _rho) {
                _p = p0;
            }
        }
        p = _p;
        return true;
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
        NT T, Lmax = _Len, max_dist, rad = 0.0, radius = P.InnerBall().second;

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
            if (max_dist > Lmax) {
                Lmax = max_dist;
            }
            pointset.push_back(p);
        }

        if (Lmax > _Len) {
            if (P.dimension() <= 2500) {
                update_delta(Lmax);
            }
        }
        pointset.clear();
    }


    inline void update_delta(NT L)
    {
        _Len = L;
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
        NT T;
        _lambdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());
        _v = GetDirection<Point>::apply(n, rng, false);
        
        do {
            _p = p;
            T = rng.sample_urdist() * _Len;
            int it = 0;

            std::pair<NT, int> pbpair
                    = P.quadratic_positive_intersect(_p, _v, _Ac, _Temp, _lambdas, _Av, _facet_prev);
            if (T <= pbpair.first || pbpair.second < 0) {
                _p += ((T * T) / (-2.0*_Temp)) *_c + (T * _v);
                _lambda_prev = T;
                return;
            }
            _lambda_prev = pbpair.first;
            _p += ((_lambda_prev * _lambda_prev) / (-2.0*_Temp)) *_c + (_lambda_prev * _v);
            _v += (-_lambda_prev/_Temp) * _c;
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, pbpair.second);

            while (it <= _rho)
            {
                std::pair<NT, int> pbpair
                        = P.quadratic_positive_intersect(_p, _v, _Ac, _Temp, _lambdas, _Av, _lambda_prev, _facet_prev);
                if (T <= pbpair.first || pbpair.second < 0) {
                    _p += ((T * T) / (-2.0*_Temp)) *_c + (T * _v);
                    _lambda_prev = T;
                    break;
                }else if (it == _rho) {
                    _lambda_prev = rng.sample_urdist() * pbpair.first;
                    _p += ((_lambda_prev * _lambda_prev) / (-2.0*_Temp)) *_c + (_lambda_prev * _v);
                    break;
                }
                _lambda_prev = pbpair.first;
                _p += ((_lambda_prev * _lambda_prev) / (-2.0*_Temp)) *_c + (_lambda_prev * _v);
                _v += (-_lambda_prev/_Temp) * _c;
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, pbpair.second);
                it++;
            }
        } while (P.is_in(_p, IN_INSIDE_BODY_TOLLERANCE) == 0);
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


    NT _Len;
    VT _Ac;
    Point _p;
    Point _v;
    Point _c;
    NT _Temp;
    NT _lambda_prev;
    int _facet_prev;
    unsigned int _rho;
    typename Point::Coeff _lambdas;
    typename Point::Coeff _Av;
};

};


#endif // RANDOM_WALKS_EXPONENTIAL_HMC_WALK_HPP

