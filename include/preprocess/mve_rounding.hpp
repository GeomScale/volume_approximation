// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef MVE_ROUNDING_HPP
#define MVE_ROUNDING_HPP

#include "mve_computation.hpp"

template <
        typename MT,
        typename VT,
        typename Polytope,
        typename Point,
        typename NT
>
std::pair< std::pair< std::pair<MT, VT>, std::pair<MT, VT> >, NT > mve_rounding(Polytope &P, 
                                                std::pair<Point, NT> InnerBall,
                                                MT &N, VT &N_shift)
{
    std::pair<std::pair<MT, VT>, bool> iter_res;
    iter_res.second = false;

    VT x0 = InnerBall.first.getCoefficients();
    MT E, L;
    unsigned int maxiter = 150, iter = 1, d = P.dimension();

    NT R = 100.0, r = 1.0, tol = std::pow(10, -6.0), reg = std::pow(10, -3.0), round_val = 1.0;

    MT T = MT::Identity(d,d);
    VT shift = VT::Zero(d);

    while (true)
    {
        iter_res = mve_computation(P.get_mat(), P.get_vec(), x0, maxiter, tol, reg);
        E = iter_res.first.first;
        E = (E + E.transpose())/2.0;
        E = E + MT::Identity(d,d)*std::pow(10,-8.0);

        Eigen::LLT<MT> lltOfA(E); // compute the Cholesky decomposition of E
        L = lltOfA.matrixL();

        Eigen::SelfAdjointEigenSolver <MT> eigensolver(L);
        r = eigensolver.eigenvalues().minCoeff();
        R = eigensolver.eigenvalues().maxCoeff();

        //std::cout<<"R/r = "<<R/r<<", iter_res.second = "<<iter_res.second<<", reg = "<<reg<<std::endl;
        if(((R <= 6.0 * r && iter_res.second)  || iter >= 20) && iter>3){
            //std::cout<<"R/r = "<<R/r<<", break"<<std::endl;
            break;
        }

        P.shift(iter_res.first.second);
        N_shift = N_shift + N*iter_res.first.second;
        N = N * L;
        shift = shift + T * iter_res.first.second;
        T = T * L;
        round_val *= L.transpose().determinant();
        P.linear_transformIt(L);

        reg = std::max(reg / 10.0, std::pow(10, -10.0));
        P.normalize();
        x0 = VT::Zero(d);

        iter++;
    }

    std::pair< std::pair< std::pair<MT, VT>, std::pair<MT, VT> >, NT > result;

    result.first.first.first = T;
    result.first.first.second = shift;

    result.first.second.first = N;
    result.first.second.second = N_shift;

    result.second = round_val;

    return result;
}

template <
        typename MT,
        typename VT,
        typename Polytope,
        typename Point,
        typename NT
>
std::pair< std::pair<MT, VT>, NT > mve_rounding(Polytope &P, std::pair<Point, NT> InnerBall)
{
    unsigned int d = P.dimension();
    MT N = MT::Identity(d,d);
    VT shift = VT::Zero(d);
    std::pair< std::pair< std::pair<MT, VT>, std::pair<MT, VT> >, NT > result = mve_rounding(P, InnerBall, N, shift);

    std::pair< std::pair<MT, VT>, NT > res;
    res.first.first = result.first.first.first;
    res.first.second = result.first.first.second;
    res.second = result.second;

    return res;
}

#endif
