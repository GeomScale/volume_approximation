// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef CONVEX_BODIES_GEN_H
#define CONVEX_BODIES_GEN_H

#include <exception>

#include "convex_bodies/convex_body.h"

#ifndef isnan
  using std::isnan;
#endif

template <class ConvexBody>
ConvexBody generate_unit_ball(unsigned int dim) {

    typedef typename ConvexBody::MT    MT;
    typedef typename ConvexBody::VT    VT;
    typedef typename ConvexBody::NT    NT;
    typedef typename ConvexBody::PointType Point;
    typedef std::function<NT(const Point&)> func;
    typedef std::function<Point(const Point&)> grad;

    func unit_ball_func = [](const Point &x) {
        return x.dot(x) - 1;
    };

    grad unit_ball_grad = [](const Point &x) {
        return 2 * x;
    };

    std::vector<func> unit_ball_funcs{unit_ball_func};
    std::vector<grad> unit_ball_grads{unit_ball_grad};

    return ConvexBody(unit_ball_funcs, unit_ball_grads, dim);
}

template <class ConvexBody>
ConvexBody generate_unit_ball_intersect_hyperplane(unsigned int dim) {

    typedef typename ConvexBody::MT    MT;
    typedef typename ConvexBody::VT    VT;
    typedef typename ConvexBody::NT    NT;
    typedef typename ConvexBody::PointType Point;
    typedef std::function<NT(const Point&)> func;
    typedef std::function<Point(const Point&)> grad;

    func unit_ball_func = [](const Point &x) {
        return x.dot(x) - 1;
    };

    func hyperplane_func = [](const Point &x) {
        return x[0] - 0.5;
    };

    grad unit_ball_grad = [](const Point &x) {
        return 2 * x;
    };

    grad hyperplane_grad = [](const Point &x) {
        Point v = Point(x.dimension());
        v.set_coord(0, 1);
        return v;
    };

    std::vector<func> unit_ball_funcs{unit_ball_func, hyperplane_func};
    std::vector<grad> unit_ball_grads{unit_ball_grad, hyperplane_grad};

    return ConvexBody(unit_ball_funcs, unit_ball_grads, dim);
}

template <class ConvexBody>
ConvexBody generate_unit_ball_intersect_logsumexp(unsigned int dim) {

    typedef typename ConvexBody::MT    MT;
    typedef typename ConvexBody::VT    VT;
    typedef typename ConvexBody::NT    NT;
    typedef typename ConvexBody::PointType Point;
    typedef std::function<NT(const Point&)> func;
    typedef std::function<Point(const Point&)> grad;

    func unit_ball_func = [](const Point &x) {
        return x.dot(x) - pow(x.dimension(), 2);
    };

    func logsumexp_func = [](const Point &x) {
        typedef typename Point::FT NT;
        NT s = 0;
        for (unsigned int i = 0; i < x.dimension(); i++) {
            s += exp(x[i]);
        }
        return log(s);
    };

    grad unit_ball_grad = [](const Point &x) {
        return 2 * x;
    };

    grad logsumexp_grad = [](const Point &x) {
        Point z(x.dimension());
        for (unsigned int i = 0; i < x.dimension(); i++) {
            z.set_coord(i, exp(x[i]));
        }
        return (1 / z.sum()) * z;
    };

    std::vector<func> unit_ball_funcs{unit_ball_func, logsumexp_func};
    std::vector<grad> unit_ball_grads{unit_ball_grad, logsumexp_grad};

    return ConvexBody(unit_ball_funcs, unit_ball_grads, dim);
}

#endif
