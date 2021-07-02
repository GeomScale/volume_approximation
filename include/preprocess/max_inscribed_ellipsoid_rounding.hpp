// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef MAX_ELLIPSOID_ROUNDING_HPP
#define MAX_ELLIPSOID_ROUNDING_HPP

#include "max_inscribed_ellipsoid.hpp"

template 
<
    typename MT,
    typename VT,
    typename NT,
    typename Polytope,
    typename Point   
>
std::tuple<MT, VT, NT> max_inscribed_ellipsoid_rounding(Polytope &P, 
                                                        Point const& InnerPoint)
{
    std::pair<std::pair<MT, VT>, bool> iter_res;
    iter_res.second = false;

    VT x0 = InnerPoint.getCoefficients();
    MT E, L;
    unsigned int maxiter = 150, iter = 1, d = P.dimension();

    NT R = 100.0, r = 1.0, tol = std::pow(10, -6.0), reg = std::pow(10, -4.0), round_val = 1.0;

    MT T = MT::Identity(d, d);
    VT shift = VT::Zero(d);

    while (true)
    {
        // compute the largest inscribed ellipsoid in P centered at x0
        iter_res = max_inscribed_ellipsoid<MT>(P.get_mat(), P.get_vec(), x0, maxiter, tol, reg);
        E = iter_res.first.first;
        E = (E + E.transpose()) / 2.0;
        E = E + MT::Identity(d, d)*std::pow(10, -8.0); //normalize E

        Eigen::LLT<MT> lltOfA(E); // compute the Cholesky decomposition of E
        L = lltOfA.matrixL();

        Eigen::SelfAdjointEigenSolver <MT> eigensolver(L);
        r = eigensolver.eigenvalues().minCoeff();
        R = eigensolver.eigenvalues().maxCoeff();

        // check the roundness of the polytope
        if(((std::abs(R / r) <= 2.3 && iter_res.second) || iter >= 20) && iter>2){
            break;
        }

        // shift polytope and apply the linear transformation on P
        P.shift(iter_res.first.second);
        shift += T * iter_res.first.second;
        T = T * L;
        round_val *= L.transpose().determinant();
        P.linear_transformIt(L);

        reg = std::max(reg / 10.0, std::pow(10, -10.0));
        P.normalize();
        x0 = VT::Zero(d);

        iter++;
    }

    std::tuple<MT, VT, NT> result = std::make_tuple(T, shift, std::abs(round_val));
    return result;
}

#endif 
