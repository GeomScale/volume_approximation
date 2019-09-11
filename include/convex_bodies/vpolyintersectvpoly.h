// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VPOLYINTERSECTVPOLY_H
#define VPOLYINTERSECTVPOLY_H

#include <iostream>
#include <iterator>
#include <vector>

template <class VPolytope>
class IntersectionOfVpoly {
public:
    typedef typename VPolytope::NT NT;
    typedef typename VPolytope::PolytopePoint PolytopePoint;
    typedef PolytopePoint Point;
    typedef typename VPolytope::MT MT;
    typedef typename VPolytope::VT VT;
    typedef typename VPolytope::rngtype RNGType;
    std::vector<Point> vecV;
    NT rad;
    VPolytope P1;
    VPolytope P2;

    IntersectionOfVpoly() {}

    IntersectionOfVpoly(VPolytope &P, VPolytope &Q) : P1(P), P2(Q) {};

    VPolytope first() { return P1; }
    VPolytope second() { return P2; }

    int is_in(Point p){
        if(P1.is_in(p)==-1)
            return P2.is_in(p);
        return 0;
    }

    void init(VPolytope &P, VPolytope &Q) {
        P1 = P;
        P2 = Q;
    }

    int num_of_hyperplanes(){
        return 0;
    }

    unsigned int dimension() {
        return P1.dimension();
    }

    int num_of_vertices() {
        return P1.num_of_vertices() + P2.num_of_vertices();
    }

    unsigned int upper_bound_of_hyperplanes() {
        return P1.upper_bound_of_hyperplanes() + P2.upper_bound_of_hyperplanes() ;
        //return 4;
    }

    std::vector<Point> get_vertices() {
        return vecV;
    }

    NT getRad() {
        return rad;
    }

    MT get_mat() {
        return P1.get_mat();
    }

    VT get_vec() {
        return P1.get_vec();
    }

    MT get_mat1() {
        return P1.get_mat();
    }

    MT get_mat2() {
        return P2.get_mat();
    }

    void set_mat(MT A) {};

    void print() {
        //std::cout<<"First polytope:\n";
        P1.print();
        //std::cout<<"\n";
        //std::cout<<"Second polytope:\n";
        P2.print();
    }

    std::pair<Point,NT> getInnerPoint_rad(bool &empty) {

        unsigned int num = 0;
        unsigned int d = P1.dimension();
        MT V1 = P1.get_mat();
        MT V2 = P2.get_mat();
        Point p(d);//, direction;
        int k1 = V1.rows();
        int k2 = V2.rows();
        int k = k1 + k2;
        Point direction(k);
        std::pair<Point, NT> cheball;
        std::vector<Point> vertices;
        typename std::vector<Point>::iterator rvert;
        bool same, done = false;

        while(true) {

            while(num<d+1){

                direction = get_direction<RNGType, Point, NT>(k);
                p = PointInIntersection<VT>(V1, V2, direction, empty);

                if (empty) {
                    return cheball;
                }

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

            cheball = P1.get_center_radius_inscribed_simplex(vertices.begin(), vertices.end(), done);
            if (done) {
                vecV = vertices;
                rad = cheball.second;
                return cheball;
            }
            vertices.clear();

            num = 0;

        }



    }



    std::pair<Point,NT> ComputeInnerBall() {

        std::pair<Point,NT> res;
        return res;

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
    std::pair<NT,NT> line_intersect(Point r,
                                    Point v) {

        std::pair <NT, NT> P1pair = P1.line_intersect(r, v);
        std::pair <NT, NT> P2pair = P2.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(P1pair.first, P2pair.first),
                                 std::max(P1pair.second, P2pair.second));

    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {
        std::pair <NT, NT> P1pair = P1.line_intersect_coord(r, rand_coord, lamdas);
        std::pair <NT, NT> P2pair = P2.line_intersect_coord(r, rand_coord, lamdas);
        return std::pair<NT, NT>(std::min(P1pair.first, P2pair.first),
                                 std::max(P1pair.second, P2pair.second));
    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          unsigned int rand_coord,
                                          unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas) {
        std::pair <NT, NT> P1pair = P1.line_intersect_coord(r, r_prev, rand_coord, rand_coord_prev, lamdas);
        std::pair <NT, NT> P2pair = P2.line_intersect_coord(r, r_prev, rand_coord, rand_coord_prev, lamdas);
        return std::pair<NT, NT>(std::min(P1pair.first, P2pair.first),
                                 std::max(P1pair.second, P2pair.second));
    }


    // shift polytope by a point c
    void shift(VT c) {
        P1.shift(c);
        P2.shift(c);
    }


    // apply linear transformation, of square matrix T, to the V-Polytope
    void linear_transformIt(MT T) {
        P1.linear_transformIt(T);
        P2.linear_transformIt(T);
    }

    std::vector<NT> get_dists(NT radius) {
        std::vector <NT> res(upper_bound_of_hyperplanes(), radius);
        return res;
    }

    template <class PointList>
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

};


#endif
