// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef SVD_ROUNDING_SINGLE_STEP_HPP
#define SVD_ROUNDING_SINGLE_STEP_HPP


template
<
    typename WalkTypePolicy,
    typename MT,
    typename VT,
    typename Polytope,
    typename Point,
    typename NT,
    typename svd_params,
    typename RandomNumberGenerator
>

void svd_rounding_single_step(Polytope &P,
                              std::pair<Point,NT> &InnerBall,
                              const unsigned int &walk_length,
                              svd_params &parameters,
                              RandomNumberGenerator &rng)
{

    int n = P.dimension(); 

    MT round_mat, r_inv;
    VT shift(n), s(n);

    Point p(n);

    NT num_its, prev_max_s = std::numeric_limits<NT>::max(),
       s_cutoff, p_cutoff;
    MT V(n,n), S(n,n);

    parameters.round_it = 1;
    parameters.max_s = std::numeric_limits<NT>::max();
    s_cutoff = 2.3;
    p_cutoff = 10.0;
    num_its = 20;

    p = InnerBall.first;
    svd_on_sample<WalkTypePolicy>(P, p, parameters.num_rounding_steps, V, s,
                                  shift, walk_length, rng);

    parameters.max_s = s.maxCoeff();

    if (parameters.max_s <= p_cutoff && parameters.max_s > s_cutoff) {
        if (parameters.last_round_under_p) {
            parameters.num_rounding_steps = parameters.num_rounding_steps * 2;
            p = InnerBall.first;
            svd_on_sample<WalkTypePolicy>(P, p, parameters.num_rounding_steps, V, s,
                                          shift, walk_length, rng);
            parameters.max_s = s.maxCoeff();
        } else {
            parameters.last_round_under_p = true;
        }
    } else {
        parameters.last_round_under_p = false;
    }
    S = s.asDiagonal();
    round_mat = V * S;
    r_inv = VT::Ones(n).cwiseProduct(s.cwiseInverse()).asDiagonal() * V.transpose();

    if (parameters.round_it != 1 && parameters.max_s >= NT(4) * parameters.prev_max_s) {
        parameters.fail = true;
        return;
    }

    parameters.round_it++;
    parameters.prev_max_s = parameters.max_s;

    P.shift(shift);
    P.linear_transformIt(round_mat);
    P.normalize();
    parameters.T_shift += parameters.T * shift;
    parameters.T = parameters.T * round_mat;

    if (parameters.max_s <= s_cutoff || parameters.round_it > num_its) {
        parameters.converged = true;
    }
}


#endif

