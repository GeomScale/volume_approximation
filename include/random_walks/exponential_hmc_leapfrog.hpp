// VolEsti (volume computation and sampling library)

// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_HAMILTONIAN_HMC_LEAPFROG_HPP
#define RANDOM_WALKS_HAMILTONIAN_HMC_LEAPFROG_HPP

#include "sampling/sphere.hpp"


// Hamiltonian Monte Carlo with leapfrog for sampling from the exponential distribution

struct HMCLeafrogExponential
{
    HMCLeafrogExponential(unsigned int steps)
            :   param(steps, true)
    {}

    HMCLeafrogExponential()
            :   param(0, false)
    {}

    struct parameters
    {
        parameters(unsigned int steps, bool set)
                :   nsteps(steps), set_steps(set)
        {}
        unsigned int nsteps;
        bool set_steps;
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
        Walk(GenericPolytope const& P, Point const& p, Point const& c, NT const& T, 
             NT const& eta, RandomNumberGenerator &rng)
        {
            _update_parameters = update_parameters();
            NT L = compute_diameter<GenericPolytope>
                    ::template compute<NT>(P);
            _steps = int(L / eta);
            _eta = eta;
            _Temp = T;
            _c = c;
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            initialize(P, p, rng);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, Point const& p, Point const& c, NT const& T, 
             NT const& eta, RandomNumberGenerator &rng, parameters const& params)
        {
            _update_parameters = update_parameters();
            _eta = eta;
            _Temp = T;
            if (params.set_steps) {
                _steps = params.nsteps;
            } else {
                NT L = compute_diameter<GenericPolytope>
                        ::template compute<NT>(P);
                _steps = int(L / _eta);
            }
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            _c = c;
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
                _v = GetDirection<Point>::apply(n, rng, false);
                _p0 = p;
                _H = _c.dot(_p) + 0.5 * _v.dot(_v);
                for (auto k=0u; k<_steps; ++k)
                {
                    T = _eta;
                    _v += (_eta / (-2.0 * _Temp)) * _c;

                    it = 0;
                    std::pair<NT, int> pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, 
                                                                          _lambda_prev, _update_parameters);
                    if (T <= pbpair.first) {
                        _p += (T * _v);
                        _v += (_eta / (-2.0 * _Temp)) * _c;
                        _lambda_prev = T;
                        continue;
                    }

                    _lambda_prev = dl * pbpair.first;
                    _p += (_lambda_prev * _v);
                    T -= _lambda_prev;
                    P.compute_reflection(_v, _p, _update_parameters);
                    it++;

                    while (it < 500*n)
                    {
                        std::pair<NT, int> pbpair
                                = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, 
                                                            _AA, _update_parameters);
                        if (T <= pbpair.first) {
                            _p += (T * _v);
                            _v += (_eta / (-2.0 * _Temp)) * _c;
                            _lambda_prev = T;
                            break;
                        }
                        _lambda_prev = dl * pbpair.first;
                        _p += (_lambda_prev * _v);
                        T -= _lambda_prev;
                        P.compute_reflection(_v, _p, _update_parameters);
                        it++;
                    }
                }
                _Htilde = _c.dot(_p) + 0.5 * _v.dot(_v);

                NT log_prob = _H - _Htilde < 0 ? _H - _Htilde : 0;
                NT u_logprob = log(rng.sample_urdist());

                if (u_logprob > log_prob) {
                    _p = _p0;
                    _Av = _Av_prev;
                    _lambda_prev = _lambda_prev_0;
                    _lambdas = _lambdas_prev;
                } else {
                    _Av_prev = _Av;
                    _lambda_prev_0 = _lambda_prev;
                    _lambdas_prev = _lambdas;
                }
            }
            p = _p;
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
            NT T = _eta;
            _lambdas.setZero(P.num_of_hyperplanes());
            _Av.setZero(P.num_of_hyperplanes());
            _p = p;
            _v = GetDirection<Point>::apply(n, rng, false);
            _v += (_eta / (-2.0 * _Temp)) * _c;

            int it = 0;

            std::pair<NT, int> pbpair
                    = P.line_first_positive_intersect(_p, _v, _lambdas, _Av, _update_parameters);
            if (T <= pbpair.first) {
                _p += (T * _v);
                _lambda_prev = T;

                _Av_prev = _Av;
                _lambda_prev_0 = _lambda_prev;
                _lambdas_prev = _lambdas;
                return;
            }
            _lambda_prev = dl * pbpair.first;
            _p += (_lambda_prev * _v);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, _update_parameters);

            while (it <= 500*n)
            {
                std::pair<NT, int> pbpair
                        = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_parameters);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    break;
                } else if (it == 500*n) {
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
            _Av_prev = _Av;
            _lambda_prev_0 = _lambda_prev;
            _lambdas_prev = _lambdas;
        }


        Point _p, _v, _c, _p0;
        NT _lambda_prev, _lambda_prev_0, _eta, _H, _Htilde, _Temp;
        unsigned int _steps;
        MT _AA;
        update_parameters _update_parameters;
        typename Point::Coeff _lambdas, _lambdas_prev;
        typename Point::Coeff _Av, _Av_prev;
    };

};


#endif



