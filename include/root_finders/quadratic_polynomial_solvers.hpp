// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef QUADRATIC_POLYNOMIAL_SOLVERS_H
#define QUADRATIC_POLYNOMIAL_SOLVERS_H

// The function compute the roots of a quadratic polynomial equation
template <typename NT>
bool solve_quadratic_polynomial(NT const& a, NT const& b, NT const& c, NT &x1, NT &x2) 
{
    if (a == NT(0)) {
        x1 = -c / b;
        x2 = x1;
        return true;
    }

    NT Delta = b * b - 4.0 * a * c;
    if (Delta < NT(0)) 
    {
        return false;
    }

    if (b >= NT(0))
    {
        x1 = (- b - std::sqrt(Delta)) / (2.0 * a);
        x2 = (2.0 * c) / (- b - std::sqrt(Delta));
    } 
    else 
    {
        x1 = (2.0 * c) / (- b + std::sqrt(Delta));
        x2 = (- b + std::sqrt(Delta)) / (2.0 * a);
    }
    return true;
}


#endif

