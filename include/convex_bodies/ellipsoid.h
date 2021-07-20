// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis
// Copyright (c) 2021- Vaibhav Thakkar

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
//Contributed and/or modified by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef ELLIPSOIDS_H
#define ELLIPSOIDS_H

#include <iostream>
#include <Eigen/Eigen>
#include <boost/math/special_functions/gamma.hpp>


template <typename NT>
NT log_gamma_function(NT x) 
{
    if (x <= NT(100)) return std::log(tgamma(x));
    return (std::log(x - NT(1)) + log_gamma_function(x - NT(1)));
}


template <class Point, class MT>
class Ellipsoid{
private:
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    // representation is (x - c)' A (x - c) <= 1
    MT A;
    Point c;

    // MT L;   // LL' = A
    unsigned int dim;

    // eigen vectors and values
    VT eigen_values;
    VT eigen_values_inv;
    VT eigen_values_inv_sqrt;
    MT eigen_vecs;

public:

    Ellipsoid(MT& Ain) : A(Ain) {
        // Eigen::LLT<Eigen::MatrixXd> lltOfA(A); // compute the Cholesky decomposition of Ain
        // if(lltOfA.info() == Eigen::NumericalIssue) {
        //     throw std::runtime_error("Possibly non semi-positive definitie matrix!");
        // }
        // L = lltOfA.matrixL();

        Eigen::SelfAdjointEigenSolver<MT> eigensolver(A);
        if (eigensolver.info() != Eigen::Success) {
            throw std::runtime_error("Eigen solver returned error!");
        }

        eigen_values = eigensolver.eigenvalues();
        eigen_vecs = eigensolver.eigenvectors();

        eigen_values_inv = eigen_values.array().inverse().matrix();
        eigen_values_inv_sqrt = eigen_values_inv.array().sqrt().matrix();

        dim = A.rows();
        c = Point(dim);
        c.set_to_origin();
    }


    Ellipsoid(MT& Ain, Point& center) : A(Ain), c(center) {
        // Eigen::LLT<Eigen::MatrixXd> lltOfA(A); // compute the Cholesky decomposition of Ain
        // if(lltOfA.info() == Eigen::NumericalIssue) {
        //     throw std::runtime_error("Possibly non semi-positive definitie matrix!");
        // }
        // L = lltOfA.matrixL();

        Eigen::SelfAdjointEigenSolver<MT> eigensolver(A);
        if (eigensolver.info() != Eigen::Success) {
            throw std::runtime_error("Eigen solver returned error!");
        }

        eigen_values = eigensolver.eigenvalues();
        eigen_vecs = eigensolver.eigenvectors();

        eigen_values_inv = eigen_values.array().inverse().matrix();
        eigen_values_inv_sqrt = eigen_values_inv.array().sqrt().matrix();

        dim = A.rows();
    }


    // Constructor for copula ellipsoid
    Ellipsoid(std::vector<std::vector<NT> >& Ain) {
        dim = Ain.size();
        A.resize(dim, dim);
        for (unsigned int i = 0; i < Ain.size(); i++) {
            for (unsigned int j = 0; j < Ain.size(); j++) {
                A(i,j) = Ain[i][j];
            }
        }

        Eigen::SelfAdjointEigenSolver<MT> eigensolver(A);
        if (eigensolver.info() != Eigen::Success) {
            throw std::runtime_error("Eigen solver returned error!");
        }

        eigen_values = eigensolver.eigenvalues();
        eigen_vecs = eigensolver.eigenvectors();

        eigen_values_inv = eigen_values.array().inverse().matrix();
        eigen_values_inv_sqrt = eigen_values_inv.array().sqrt().matrix();

        dim = A.rows();
        c = Point(dim);
        c.set_to_origin();
    }


    VT eigenvals() const {
        return eigen_values;
    }


    VT eigenvals_inv() const {
        return eigen_values_inv;
    }


    VT eigenvals_inv_sqrt() const {
        return eigen_values_inv_sqrt;
    }


    MT eigenvecs() const {
        return eigen_vecs;
    }


    unsigned int dimensions() const {
        return dim;
    }


    void print() {
        std::cout << "Ellipse is in the form: (x-c)' A (x-c) <= 1 \n";
        std::cout << "c = \n" << c.print();
        std::cout << "A = \n" << A;
    }



    NT mat_mult(Point const& p) {
        NT res = NT(0);
        for (size_t i=0; i<dim; i++) {
            for (size_t j = i; j<dim; j++) {
                res += 2*A(i, j)*p[i]*p[j];
            }
        }

        return res;
    }

    
    NT log_volume () {
        NT ball_log_vol = (NT(dim)/NT(2) * std::log(M_PI)) - log_gamma_function(NT(dim) / NT(2) + 1);
        NT det_factor = - 0.5 * std::log( A.determinant() );

        return det_factor + ball_log_vol;
    }


    void scale(NT scale_factor) {
        assert (scale_factor > 0);

        NT scale_factor_sq = scale_factor * scale_factor;
        NT inv_scale_factor = (NT(1.0) / scale_factor);
        NT inv_scale_factor_sq = (NT(1.0) / scale_factor_sq);

        // L = mult_factor * L;
        eigen_values = inv_scale_factor_sq * eigen_values;
        eigen_values_inv = scale_factor_sq * eigen_values_inv;
        eigen_values_inv_sqrt = scale_factor * eigen_values_inv_sqrt;

        A = inv_scale_factor_sq * A; // as volume depends on square root of it's determinant
    }


    int is_in(Point const& p){
        NT val = mat_mult(p - c);
        if (val > 1) {
            return 0;
        }

        return -1;
    }


    // compute intersection point of ray starting from r and pointing to v
    std::pair<NT, NT> line_intersect(Point& r, Point& v) const {
        // constants of a quadratic equation
        NT a_q = mat_mult(v);
        NT b_q = 2 * (r - c).getCoefficients().dot(A * v.getCoefficients());
        NT c_q = mat_mult(r - c);
        
        NT D = std::pow(b_q, 2) - 4*a_q*c_q;
        return std::pair<NT, NT> ( (-b_q + std::sqrt(D))/(2*a_q) , (-b_q - std::sqrt(D))/(2*a_q) );
    }


    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    const VT &Ar,
                                    const VT &Av) const
    {
        return line_intersect(r, v);
    }


    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    const VT &Ar,
                                    const VT &Av,
                                    NT &lambda_prev) const
    {
        return line_intersect(r, v);
    }


    std::pair<NT,int> line_positive_intersect(Point const& r,
                                              Point const& v) const
    {
        return std::pair<NT,NT>(line_intersect(r, v).first, 0);
    }


    std::pair<NT,int> line_positive_intersect(Point const& r,
                                              Point const& v,
                                              const VT &Ar,
                                              const VT &Av) const
    {
        return line_positive_intersect(r, v);
    }


    std::pair<NT,int> line_positive_intersect(Point const& r,
                                              Point const& v,
                                              const VT &Ar,
                                              const VT &Av,
                                              NT &lambda_prev) const
    {
        return line_positive_intersect(r, v);
    }


    // Compute the intersection of a coordinate ray
    std::pair<NT,NT> line_intersect_coord(const Point &r, const unsigned int rand_coord) const {
        NT a_q = A(rand_coord, rand_coord);
        NT b_q = 2 * (r - c).getCoefficients().dot(A.col(rand_coord));
        NT c_q = mat_mult(r - c);
        
        NT D = std::pow(b_q, 2) - 4*a_q*c_q;
        return std::pair<NT, NT> ( (-b_q + std::sqrt(D))/(2*a_q) , (-b_q - std::sqrt(D))/(2*a_q) );
    }


    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          unsigned int const& rand_coord,
                                          const VT &lamdas) const
    {
        return line_intersect_coord(r, rand_coord);
    }


    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          Point const& r_prev,
                                          unsigned int const& rand_coord,
                                          unsigned int const& rand_coord_prev,
                                          const VT &lamdas) const
    {
        return line_intersect_coord(r, rand_coord);
    }


    void compute_reflection (Point& v, Point const& p) const
    {
        // normal vector is Ap
        Point s(A * (p-c).getCoefficients());
        s *= (1.0 / std::sqrt(s.squared_length()));
        s *= (-2.0 * v.dot(s));
        v += s;
    }
};

#endif
