// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018-19 programs.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VPOLYINTERSECTVPOLY_H
#define VPOLYINTERSECTVPOLY_H

#include <iostream>
#include <iterator>
#include <vector>
#include "sampling/sphere.hpp"


template <typename VPolytope, typename RNGType>
class IntersectionOfVpoly {
public:
    typedef typename VPolytope::NT NT;
    typedef typename VPolytope::PointType PointType;
    typedef PointType Point;
    typedef typename VPolytope::MT MT;
    typedef typename VPolytope::VT VT;
    unsigned seed;
    std::pair<Point,NT> _inner_ball;
    NT rad;
    VPolytope P1;
    VPolytope P2;

    IntersectionOfVpoly(): P1(), P2() {}

    IntersectionOfVpoly(VPolytope P, VPolytope Q) : P1(P), P2(Q) {
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    IntersectionOfVpoly(VPolytope P, VPolytope Q, unsigned _seed) : P1(P), P2(Q)  {
        seed = _seed;
    }

    VPolytope first() { return P1; }
    VPolytope second() { return P2; }

    std::pair<Point,NT> InnerBall() const
    {
        return _inner_ball;
    }

    int is_in(const Point &p) const {
        if(P1.is_in(p)==-1)
            return P2.is_in(p);
        return 0;
    }


    int num_of_hyperplanes() const {
        return 0;
    }

    unsigned int dimension() const {
        return P1.dimension();
    }

    int num_of_vertices() const {
        return P1.num_of_vertices() + P2.num_of_vertices();
    }

    unsigned int upper_bound_of_hyperplanes() const {
        return dimension() + 1;
        //return 4;
    }

    int num_of_generators() const {
        return 0;
    }

    NT getRad() const {
        return rad;
    }

    MT get_mat() const {
        return P1.get_mat();
    }

    VT get_vec() const {
        return P1.get_vec();
    }

    MT get_T() const {
        return P1.get_mat();
    }

    MT get_mat2() const {
        return P2.get_mat();
    }

    Point get_mean_of_vertices() const {
        return Point(P1.dimension());
    }


    NT get_max_vert_norm() const {
        return 0.0;
    }

    void print() {
        P1.print();
        P2.print();
    }

    bool is_feasible() {
        bool empty;
        int k = P1.get_mat().rows() + P2.get_mat().rows();
        RNGType rng(k);
        rng.set_seed(seed);
        PointInIntersection<VT>(P1.get_mat(), P2.get_mat(),
                                GetDirection<Point>::apply(k, rng), empty);
        return !empty;
    }

    std::pair<Point,NT> ComputeInnerBall() {

        unsigned int num = 0, d = P1.dimension();
        MT V1 = P1.get_mat(), V2 = P2.get_mat();
        int k1 = V1.rows(), k2 = V2.rows();
        int k = k1 + k2;
        RNGType rng(k);
        rng.set_seed(seed);
        Point direction(k), p(d);
        std::vector<Point> vertices;
        typename std::vector<Point>::iterator rvert;
        bool same;

        while(num<d+1){

            direction = GetDirection<Point>::apply(k, rng);
            p = PointInIntersection<VT>(V1, V2, direction, same);

            same = false;
            rvert = vertices.begin();
            for ( ;  rvert!=vertices.end(); ++rvert) {
                if (p==(*rvert)) {
                    same = true;
                    break;
                }
            }
            if (same) continue;
            vertices.push_back(p);
            num++;

        }

        _inner_ball = P1.get_center_radius_inscribed_simplex(vertices.begin(), vertices.end());
        return _inner_ball;

    }

    void set_InnerBall(std::pair<Point,NT> const& innerball) const
    {
        _inner_ball = innerball;
    }

/*
        unsigned int num_of_v = 0;
        unsigned int d = dimension();
        MT V(0, d);
        MT V1 = P1.get_mat();
        MT V2 = P2.get_mat();
        VT itervec(d);
        std::vector<NT> temp_p(d, 0.0);
        typename std::vector<NT>::iterator pit;
        int j;
        for (int i = 0; i < V1.rows(); ++i) {
            pit = temp_p.begin();
            j = 0;
            for (; pit!=temp_p.end(); ++pit, ++j) {
                *pit = V1(i,j);
                itervec(j) = V1(i,j);
            }
            Point p(d, temp_p.begin(), temp_p.end());
            if (P2.is_in(p) == -1) {
                V.conservativeResize(V.rows() + 1, V.cols());
                V.row(V.rows()-1) = itervec;
                num_of_v++;
            }
        }

        for (int i = 0; i < V2.rows(); ++i) {
            pit = temp_p.begin();
            j = 0;
            for (; pit!=temp_p.end(); ++pit, ++j) {
                *pit = V2(i,j);
                itervec(j) = V2(i,j);
            }
            Point p(d, temp_p.begin(), temp_p.end());
            if (P1.is_in(p) == -1) {
                V.conservativeResize(V.rows() + 1, V.cols());
                V.row(V.rows()-1) = itervec;
                num_of_v++;
            }
        }
        if (num_of_v <= d) {
            std::cout<<"no simplex"<<std::endl;
            std::pair<Point,NT> res;
            res.second = -1.0;
            return res;
        }

        VPolytope Q;
        Q.init(d, V, itervec);
        return Q.ComputeInnerBall();
    }*/

    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(const Point &r, const Point &v) const {

        std::pair <NT, NT> P1pair = P1.line_intersect(r, v);
        std::pair <NT, NT> P2pair = P2.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(P1pair.first, P2pair.first),
                                 std::max(P1pair.second, P2pair.second));

    }

    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(const Point &r, const Point &v, const VT &Ar,
            const VT &Av) const {
        return line_intersect(r, v);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(const Point &r, const Point &v, const VT &Ar,
                                    const VT &Av, const NT &lambda) const {
        return line_intersect(r, v);
    }

    std::pair<NT, int> line_positive_intersect(const Point &r, const Point &v) const {

        std::pair<NT, int> P1pair = P1.line_positive_intersect(r, v);
        std::pair<NT, int> P2pair = P2.line_positive_intersect(r, v);

        if(P1pair.first < P2pair.first) {
            return std::pair<NT, int>(P1pair.first, 1);
        }
        return std::pair<NT, int>(P2pair.first, 2);
    }

    std::pair<NT, int> line_positive_intersect(const Point &r, const Point &v, const VT &Ar,
            const VT &Av) const {
        return line_positive_intersect(r, v);
    }


    std::pair<NT, int> line_positive_intersect(const Point &r, const Point &v, const VT &Ar,
                                               const VT &Av, const NT &lambda_prev) const {
        return line_positive_intersect(r, v);//, Ar, Av);
    }

    //------------------------------accelarated billiard------------------------------//
    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(Point const& r,
                                                     Point const& v,
                                                     VT& Ar,
                                                     VT& Av,
                                                     update_parameters& params) const
    {
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               MT const& AA,
                                               update_parameters& params) const
    {
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               update_parameters& params) const
    {
        return line_positive_intersect(r, v);
    }
    //------------------------------------------------------------------------------//


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(const Point &r,
                                          const unsigned int &rand_coord,
                                          const VT &lamdas) const {
        std::pair <NT, NT> P1pair = P1.line_intersect_coord(r, rand_coord, lamdas);
        std::pair <NT, NT> P2pair = P2.line_intersect_coord(r, rand_coord, lamdas);
        return std::pair<NT, NT>(std::min(P1pair.first, P2pair.first),
                                 std::max(P1pair.second, P2pair.second));
    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(const Point &r,
                                          const Point &r_prev,
                                          const unsigned int &rand_coord,
                                          const unsigned int &rand_coord_prev,
                                          const VT &lamdas) const {
        return line_intersect_coord(r, rand_coord, lamdas);
    }


    // shift polytope by a point c
    void shift(const VT &c) {
        P1.shift(c);
        P2.shift(c);
    }


    // apply linear transformation, of square matrix T, to the V-Polytope
    void linear_transformIt(const MT &T) {
        P1.linear_transformIt(T);
        P2.linear_transformIt(T);
    }

    std::vector<NT> get_dists(const NT &radius) const {
        std::vector <NT> res(upper_bound_of_hyperplanes(), radius);
        return res;
    }

    template <typename PointList>
    bool get_points_for_rounding (PointList &randPoints) {
        if (num_of_vertices()>40*dimension()) {
            return false;
        }
        if(!P1.get_points_for_rounding(randPoints)) {
            return false;
        }
        if(!P2.get_points_for_rounding(randPoints)) {
            return false;
        }

        return true;
    }

    void normalize() {}

    void resetFlags() {}

    void update_position_internal(NT&){}
    
    void compute_reflection (Point &v, const Point &p, const int &facet) const {

        if (facet == 1) {
            P1.compute_reflection (v, p, facet);
        } else {
            P1.compute_reflection (v, p, facet);
        }

    }

    template <typename update_parameters>
    void compute_reflection (Point &v, const Point &p, update_parameters const& params) const {

        if (params.facet_prev == 1) {
            P1.compute_reflection (v, p, params);
        } else {
            P2.compute_reflection (v, p, params);
        }

    }

};


#endif
