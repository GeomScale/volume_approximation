// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis


// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GAUSSIAN_EXACT_HMC_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_EXACT_HMC_WALK_HPP

#include "sampling/sphere.hpp"



// Exact HMC for sampling from the spherical Gaussian distribution

struct GaussianHamiltonianMonteCarloExactWalk
{
    GaussianHamiltonianMonteCarloExactWalk(double L, unsigned int _rho)
            :   param(L, true, _rho, true)
    {}
    
    GaussianHamiltonianMonteCarloExactWalk(double L)
            :   param(L, true, 0, false)
    {}

    GaussianHamiltonianMonteCarloExactWalk()
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
    typedef typename Polytope::VT VT;

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p, NT const& a_i, RandomNumberGenerator &rng)
    {
        _Len = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
        _omega = std::sqrt(NT(2) * a_i);
        _rho = 100 * P.dimension(); // upper bound for the number of reflections (experimental)
        initialize(P, p, a_i, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p, NT const& a_i, RandomNumberGenerator &rng,
         parameters const& params)
    {
        _Len = params.set_L ? params.m_L
                          : compute_diameter<GenericPolytope>
                            ::template compute<NT>(P);
        _omega = std::sqrt(NT(2) * a_i);
        _rho = 100 * P.dimension(); // upper bound for the number of reflections (experimental)
        initialize(P, p, a_i, rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline void apply(GenericPolytope const& P,
                      Point& p,
                      NT const& a_i,
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT T;

        for (auto j=0u; j<walk_length; ++j)
        {
            T = rng.sample_urdist() * _Len;
            _v = GetDirection<Point>::apply(n, rng, false);
            Point p0 = _p;
            int it = 0;
            while (it < _rho)
            {
                auto pbpair = P.trigonometric_positive_intersect(_p, _v, _omega, _facet_prev);
                if (T <= pbpair.first) {
                    update_position(_p, _v, T, _omega);
                    break;
                }
                _lambda_prev = pbpair.first;
                T -= _lambda_prev;
                update_position(_p, _v, _lambda_prev, _omega);
                P.compute_reflection(_v, _p, pbpair.second);
                it++;
            }
            if (it == _rho){
                _p = p0;
            }
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
                           NT const& a_i,
                           RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        _facet_prev = -1;
        _p = p;
        _v = GetDirection<Point>::apply(n, rng, false);

        NT T = rng.sample_urdist() * _Len;
        int it = 0;

        while (it <= _rho)
        {
            auto pbpair
                    = P.trigonometric_positive_intersect(_p, _v, _omega, _facet_prev);
            if (T <= pbpair.first) {
                update_position(_p, _v, T, _omega);
                break;
            }else if (it == _rho) {
                _lambda_prev = rng.sample_urdist() * pbpair.first;
                update_position(_p, _v, _lambda_prev, _omega);
                break;
            }
            _lambda_prev = pbpair.first;
            update_position(_p, _v, _lambda_prev, _omega);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, pbpair.second);
            it++;
        }
    }

    inline void update_position(Point &p, Point &v, NT const& T, NT const& omega)
    {
        NT C, Phi;
        for (size_t i = 0; i < p.dimension(); i++)
        {
            C = std::sqrt(p[i] * p[i] + (v[i] * v[i]) / (omega * omega));
            Phi = std::atan((-v[i]) / (p[i] * omega));
            if (v[i] < 0.0 && Phi < 0.0) {
                Phi += M_PI;
            } else if (v[i] > 0.0 && Phi > 0.0) {
                Phi -= M_PI;
            }
            p.set_coord(i, C * std::cos(omega * T + Phi));
            v.set_coord(i, -C * omega * std::sin(omega * T + Phi));
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

    int _facet_prev;
    unsigned int _rho;
    NT _Len;
    Point _p;
    Point _v;
    NT _omega;
    NT _lambda_prev;
};

};


#endif // RANDOM_WALKS_GAUSSIAN_HMC_WALK_HPP

