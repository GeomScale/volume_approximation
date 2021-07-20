// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2021- Vaibhav Thakkar

//Contributed and/or modified by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.


#ifndef SAMPLERS_ELLIPSOID_HPP
#define SAMPLERS_ELLIPSOID_HPP

#include "sphere.hpp"
#include "Eigen/Eigen"


template <typename Point>
struct GetPointInDellipsoid
{
    template <typename NT, typename VT, typename MT, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              VT const& eigenvals_inv_sqrt, // sqrt of inverse of eigenvals of matrix A in (x'Ax <= 1)
                              MT const& EigenVecs,          // eigenvecs of matrix A in (x'Ax <= 1)
                              RandomNumberGenerator &rng)
    {
        // Generate a point inside a sphere of radius 1.0
        Point p = GetPointInDsphere<Point>::apply(dim, NT(1.0), rng);

        // scale points to the ellipsoid using the eigenvalues
        VT scaled_vec = p.getCoefficients().cwiseProduct(eigenvals_inv_sqrt);

        // rotate with the eigenvectors
        return Point(EigenVecs * scaled_vec);
    }
};



#endif // ELLIPSOID_HPP
