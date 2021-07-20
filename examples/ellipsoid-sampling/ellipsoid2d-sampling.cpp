#include <iostream>
#include <boost/random.hpp>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "ellipsoid.h"
#include "sampling/ellipsoid.hpp"
#include "random_walks/random_walks.hpp"

typedef double NT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

int main(int argc, char const *argv[]) {
    unsigned int dim = 2;
    MT L(2, 2);
    L << 0.5, 0,
         1.5, 1.0;
    MT A = L * L.transpose();
    Point c(dim, {2, 2});

    Ellipsoid<Point, MT> ell(A, c);
    
    int num_points = 1000;
    Point p(dim);
    RNGType rng(dim);
    VT eigenvals_inv_sqrt = ell.eigenvals_inv_sqrt();
    MT eigenvecs = ell.eigenvecs();

    for (int i=0; i<num_points; ++i) {
        p = GetPointInDellipsoid<Point>::apply<NT>(dim, eigenvals_inv_sqrt, eigenvecs, rng);
        p += c;
        p.print();
    }
    return 0;
}

