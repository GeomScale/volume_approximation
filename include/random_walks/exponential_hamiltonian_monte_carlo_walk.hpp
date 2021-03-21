// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis


// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_EXPONENTIAL_HMC_WALK_HPP
#define RANDOM_WALKS_EXPONENTIAL_HMC_WALK_HPP

#include "sampling/sphere.hpp"


// Exact HMC for sampling from the Exponential distribution restricted to a convex polytope

struct HMCExponentialWalk
{
    HMCExponentialWalk(double L)
            :   param(L, true)
    {}

    HMCExponentialWalk()
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
        Point p0;

        for (auto j=0u; j<walk_length; ++j)
        {
            do {
                failures++;
                if (failures == 1000) {
                    return false;
                }
                T = -std::log(rng.sample_urdist()) * _Len;
                _v = GetDirection<Point>::apply(n, rng, false);
                p0 = _p;
                it = 0;
                while (it < 200*n)
                {
                    auto pbpair = P.quadratic_positive_intersect(_p, _v, _Ac, _Temp, _lambdas,
                                                             _Av, _lambda_prev, _facet_prev);
                    if (T <= pbpair.first || pbpair.second < 0) {
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
                
            } while (P.is_in(_p, _tol) == 0);
            if (it == 200*n){
                _p = p0;
            }
        }
        p = _p;
        return true;
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
        _p = p;
        _v = GetDirection<Point>::apply(n, rng, false);
        
        do {
            T = -std::log(rng.sample_urdist()) * _Len;
            Point p0 = _p;
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

            while (it <= 100*n)
            {
                std::pair<NT, int> pbpair
                        = P.quadratic_positive_intersect(_p, _v, _Ac, _Temp, _lambdas, _Av, _lambda_prev, _facet_prev);
                if (T <= pbpair.first || pbpair.second < 0) {
                    _p += ((T * T) / (-2.0*_Temp)) *_c + (T * _v);
                    _lambda_prev = T;
                    break;
                }else if (it == 100*n) {
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
        } while (P.is_in(_p, _tol) == 0);
    }

    NT _Len;
    VT _Ac;
    Point _p;
    Point _v;
    Point _c;
    NT _Temp;
    const double _tol = 1e-10;
    NT _lambda_prev;
    int _facet_prev;
    typename Point::Coeff _lambdas;
    typename Point::Coeff _Av;
};

};


#endif // RANDOM_WALKS_EXPONENTIAL_HMC_WALK_HPP

