// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis


// Licensed under GNU LGPL.3, see LICENCE file


#ifndef QUADRATIC_POLYNOMIAL_SOLVERS_H
#define QUADRATIC_POLYNOMIAL_SOLVERS_H

// The function returns the sigh of the input number
template <typename T> int sgn(T val) 
{
    return (T(0) < val) - (val < T(0));
}

// The function compute the roots of a quadratic polynomial equation
template <typename NT>
void solve_quadratic_polynomial(NT a, NT b, NT c, NT &x1, NT &x2, bool &real) {

    NT Delta = b * b - 4.0 * a * c;
    if (Delta < NT(0)) {
        real = false;
        return;
    }
    real = true;
    if (b >= NT(0)){
        x1 = (- b - std::sqrt(Delta)) / (2.0 * a);
        x2 = (2.0 * c) / (- b - std::sqrt(Delta));
    } else {
        x1 = (2.0 * c) / (- b + std::sqrt(Delta));
        x2 = (- b + std::sqrt(Delta)) / (2.0 * a);
    }
}

// The function compute the roots of a quadratic polynomial equation
// using the algorithm in Revisiting the stability of computing the roots of a quadratic polynomial
// by Nicola Mastronardi and Paul Michel Van Dooren, Electronic transactions on numerical analysis, 2014
template <typename NT>
void solve_qudratic_polynomial_stable(NT a, NT b, NT c, NT &x1, NT &x2, bool &real)
{
    real = true;

    if (c == NT(0)) {
        x1 = NT(0);
        x2 = -b/a;
        return;
    }
    b = b / a;
    c = c / a;
    
    NT alpha = sgn(b) * std::sqrt(std::abs(c));
    NT e = NT(sgn(c));
    NT beta = b / (NT(2)*alpha);

    if (e > NT(0) && beta < NT(1))
    {
        real = false;
        return;
    }

    if (e < NT(0))
    {
        x1 = beta + std::sqrt(beta*beta + NT(1));
        x2 = -NT(1)/x1;
    } else 
    {
        x1 = beta + std::sqrt((beta + NT(1)) * (beta - NT(1)));
        x2 = NT(1) / x1;
    }
    x1 = -x1 * alpha;
    x2 = -x2 * alpha;
}

#endif

