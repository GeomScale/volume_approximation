#ifndef CONVEX_BODY_H
#define CONVEX_BODY_H

#include <limits>
#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <functional>

template <typename Point>
class ConvexBody {
public:
    typedef Point                                             PointType;
    typedef typename Point::FT                                NT;
    typedef typename std::vector<NT>::iterator                viterator;
    //using RowMatrixXd = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    //typedef RowMatrixXd MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;
    typedef std::function<NT(const Point&)>                          func;
    typedef std::function<Point(const Point&)>                       grad;
private:
    unsigned int        dim; //dimension
    std::vector<func>   gs; // convex functions defining the convex body
    std::vector<grad>   grad_gs; // convex function gradients
    unsigned int m;
    NT tol = NT(1e-4);


public:

    ConvexBody() : m(0) {}

    ConvexBody(std::vector<func> gs_, std::vector<grad> grad_gs_, unsigned int dim_) :
        gs(gs_), grad_gs(grad_gs_), dim(dim_)
    {
        m = gs.size();
    }

    unsigned int dimension() {
        return dim;
    }

    // Compute positive line intersection (in [0, 1]) using Binary search
    // x: starting point
    // v: direction (ray)
    std::pair<NT, int> line_positive_intersect(Point const& x, Point const &v) const {
        NT t_min = NT(1);
        int constraint = -1;
        NT t;
        for (unsigned int i = 0; i < m; i++) {
            t = binary_search(x, v, gs[i]);
            if (t < t_min) {
                t_min = t;
                constraint = i;
            }
        }

        return std::make_pair(t_min, constraint);
    }

    NT binary_search(Point const &x, Point const &v, func const& f) const {
        NT t_min = NT(0);
        NT t_max = NT(1);
        NT t;
        NT value;

        while (t_max - t_min > tol) {
            t = (t_max + t_min) / 2;
            value = f(x + t * v);
            if (value >= -tol && value <= 0) {
                return t;
            } else if (value < -tol) {
                t_min = t;
            } else {
                t_max = t;
            }

        }

        return t;
    }


    // Computes unit normal at point p of the boundary
    Point unit_normal(Point const& p, int const& constraint) const {
        Point n = grad_gs[constraint](p);
        return (1 / n.length()) * n;
    }

    // Computes reflection of v about point p on the boundary
    void compute_reflection(Point &v, Point const& p, int const& constraint) const {
        Point n = unit_normal(p, constraint);
        v += -2 * v.dot(n) * n;
    }

    // Check if point is in K
    int is_in(Point const& p, NT tol=NT(0)) {
        for (func g : gs) {
            if (g(p) > NT(-tol)) return 0;
        }
        return -1;
    }

};

#endif
