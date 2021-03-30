// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file


#include <Rcpp.h>

class Spectrahedron {
public:
    /// A list with the matrices A0, ..., An
    Rcpp::List matrices;

    Spectrahedron() {}

    Spectrahedron(Rcpp::List _matrices) : matrices(_matrices) {}
};


RCPP_MODULE(spectrahedron){
        using namespace Rcpp ;

        //' An exposed class to represent a spectrahedron
        //'
        //' @description A spectrahedron is a convex body defined by a linear matrix inequality of the form \eqn{A_0 + x_1 A_1 + ... + x_n A_n \preceq 0}.
        //' The matrices \eqn{A_i} are symmetric \eqn{m \times m} real matrices and \eqn{\preceq 0} denoted negative semidefiniteness.
        //'
        //' @field matrices The matrices \eqn{A_0, A_1, ..., A_n}
        //'
        //' @example
        //' A0 = matrix(c(-1,0,0,0,-2,1,0,1,-2), nrow=3, ncol=3, byrow = TRUE)
        //' A1 = matrix(c(-1,0,0,0,0,1,0,1,0), nrow=3, ncol=3, byrow = TRUE)
        //' A2 = matrix(c(0,0,-1,0,0,0,-1,0,0), nrow=3, ncol=3, byrow = TRUE)
        //' lmi = list(M0, M1,M2)
        //' S = Spectrahedron$new(lmi);
        //'
        //' @export
        class_<Spectrahedron>("Spectrahedron")
        // expose the default constructor
        .constructor()
        .constructor<Rcpp::List>()
        .field( "matrices", &Spectrahedron::matrices);
}

extern SEXP _rcpp_module_boot_spectrahedron(void);
