// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018-19 programs.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VPOLYTOPE_H
#define VPOLYTOPE_H

#include <limits>
#include <iostream>
#include <Eigen/Eigen>

#include "lp_oracles/vpolyoracles.h"
#include <khach.h>


// V-Polytope class
template<typename Point>
class VPolytope {
public:
    typedef Point                                             PointType;
    typedef typename Point::FT                                NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;

private:
    unsigned int         _d;  //dimension
    MT                   V;  //matrix V. Each row contains a vertex
    VT                   b;  // vector b that contains first column of ine file
    std::pair<Point, NT> _inner_ball;

    // TODO: Why don't we use std::vector<REAL>  and std::vector<int> for these pointers?
    REAL *conv_comb, *conv_comb2, *conv_mem, *row;
    int *colno, *colno_mem;

public:
    VPolytope() {}

    VPolytope(const unsigned int &dim, const MT &_V, const VT &_b):
            _d{dim}, V{_V}, b{_b},
            conv_comb{new REAL[V.rows() + 1]},
            conv_comb2{new REAL[V.rows() + 1]},
            conv_mem{new REAL[V.rows()]},
            row{new REAL[V.rows() + 1]},
            colno{new int[V.rows() + 1]},
            colno_mem{new int[V.rows()]} 
    {
    }

    // Construct matrix V which contains the vertices row-wise
    // TODO: change rows;
    VPolytope(std::vector<std::vector<NT>> const& Pin)
    {
        _d = Pin[0][1] - 1;
        V.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (unsigned int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            unsigned int j;
            for (j = 1; j < _d + 1; j++) {
                V(i - 1, j - 1) = Pin[i][j];
            }
        }
        conv_comb = new REAL[Pin.size()];
        conv_comb2 = new REAL[Pin.size()];
        conv_mem = new REAL[V.rows()];
        row = new REAL[V.rows() + 1];
        colno = new int[V.rows() + 1];
        colno_mem = new int[V.rows()];
    }

    template <typename T>
    void copy_array(T* source, T* result, size_t count)
    {
        T* tarray;
        tarray = new T[count];
        std::copy_n(source, count, tarray);
        delete [] result;
        result = tarray;
    }

    VPolytope& operator=(const VPolytope& other)
    {
        if (this != &other) { // protect against invalid self-assignment
            _d = other._d;
            V = other.V;
            b = other.b;

            copy_array(other.conv_comb, conv_comb, V.rows() + 1);
            copy_array(other.conv_comb2, conv_comb2, V.rows() + 1);
            copy_array(other.conv_mem, conv_mem, V.rows());
            copy_array(other.row, row, V.rows() + 1);
            copy_array(other.colno, colno, V.rows() + 1);
            copy_array(other.colno_mem, colno_mem, V.rows());
        }
        return *this;
    }

    VPolytope& operator=(VPolytope&& other)
    {
        if (this != &other) { // protect against invalid self-assignment
            _d = other._d;
            V = other.V;
            b = other.b;

            conv_comb = other.conv_comb;  other.conv_comb = nullptr;
            conv_comb2 = other.conv_comb2;  other.conv_comb2 = nullptr;
            conv_mem = other.conv_mem;  other.conv_mem = nullptr;
            row = other.row; other.row = nullptr;
            colno = other.colno; colno = nullptr;
            colno_mem = other.colno_mem; colno_mem = nullptr;
        }
        return *this;
    }


    VPolytope(const VPolytope& other) :
            _d{other._d}, V{other.V}, b{other.b},
            conv_comb{new REAL[V.rows() + 1]},
            conv_comb2{new REAL[V.rows() + 1]},
            conv_mem{new REAL[V.rows()]},
            row{new REAL[V.rows() + 1]},
            colno{new int[V.rows() + 1]},
            colno_mem{new int[V.rows()]}
    {
        std::copy_n(other.conv_comb, V.rows() + 1, conv_comb);
        std::copy_n(other.conv_comb2, V.rows() + 1, conv_comb2);
        std::copy_n(other.conv_mem, V.rows(), conv_mem);
        std::copy_n(other.row, V.rows() + 1, row);
        std::copy_n(other.colno, V.rows() + 1, colno);
        std::copy_n(other.colno_mem, V.rows(), colno_mem);
    }

    VPolytope(VPolytope&& other) :
            _d{other._d}, V{other.V}, b{other.b},
            conv_comb{nullptr}, conv_comb2{nullptr}, conv_mem{nullptr}, row{nullptr},
            colno{nullptr}, colno_mem{nullptr}
    {
        conv_comb = other.conv_comb;  other.conv_comb = nullptr;
        conv_comb2 = other.conv_comb2;  other.conv_comb2 = nullptr;
        conv_mem = other.conv_mem;  other.conv_mem = nullptr;
        row = other.row; other.row = nullptr;
        colno = other.colno; colno = nullptr;
        colno_mem = other.colno_mem; colno_mem = nullptr;
    }

    ~VPolytope() {
        delete [] conv_comb;
        delete [] conv_comb2;
        delete [] colno;
        delete [] colno_mem;
        delete [] row;
        delete [] conv_mem;
    }

    std::pair<Point,NT> InnerBall() const
    {
        return _inner_ball;
    }

    void set_InnerBall(std::pair<Point,NT> const& innerball) //const
    {
        _inner_ball = innerball;
    }

    // return dimension
    unsigned int dimension() const {
        return _d;
    }


    // this function returns 0. The main sampler requests this function to set the length of lambdas vector
    int num_of_hyperplanes() const {
        return 0;
    }

    int num_of_generators() const {
        return 0;
    }

    // compute the number of facets of the cyclic polytope in dimension _d with the same number of vertices
    // this is an upper bound for the number of the facets from McMullen's Upper Bound Theorem
    unsigned int upper_bound_of_hyperplanes() const {
        return 2 * _d;
    }


    // return the number of vertices
    int num_of_vertices() const {
        return V.rows();
    }

    // return the matrix V
    MT get_mat() const {
        return V;
    }

    // return the vector b
    VT get_vec() const {
        return b;
    }

    // change the matrix V
    void set_mat(const MT &V2) {
        V = V2;
    }

    // change the vector b
    void set_vec(const VT &b2) {
        b = b2;
    }

    MT get_T() const {
        return V;
    }


    // print polytope in input format
    void print() {
        std::cout << " " << V.rows() << " " << _d << " float" << std::endl;
        for (unsigned int i = 0; i < V.rows(); i++) {
            for (unsigned int j = 0; j < _d; j++) {
                std::cout << V(i, j) << " ";
            }
            std::cout<<"\n";
        }
    }

    Point get_mean_of_vertices() {
        Point xc(_d);
        for (int i = 0; i < num_of_vertices(); ++i) {
            xc.add(V.row(i));
        }
        xc *= (1.0/NT(num_of_vertices()));

        return xc;
    }


    NT get_max_vert_norm() {
        NT rad =0.0;
        NT rad_iter;
        for (int i = 0; i < num_of_vertices(); ++i) {
            rad_iter = V.row(i).norm();
            if(rad_iter>rad)rad = rad_iter;
        }
        return rad;
    }

    void normalize() {}

    // take d+1 points as input and compute the chebychev ball of the defined simplex
    // done is true when the simplex is full dimensional and false if it is not
    std::pair<Point,NT> get_center_radius_inscribed_simplex(const typename std::vector<Point>::iterator it_beg,
                                                            const typename std::vector<Point>::iterator it_end) {

        Point p0 = *it_beg,p1,c;
        unsigned int dim = p0.dimension(), i, j;
        std::vector <NT> temp_p;
        NT radius = 0.0, gi, sum = 0.0;
        MT B(dim,dim), Bg(dim,dim), e(1,dim);
        VT row(dim), g(dim);
        std::pair <Point,NT> result;

        for (j=1; j<dim+1; j++) {
            Point* pk = &(*(it_beg + j));
            e(j - 1) = 1.0;
            B.col(j-1) = pk->getCoefficients() - p0.getCoefficients();
        }

        Bg = B;
        B = B.inverse();
        for (i=0; i<dim; i++) {
            gi = B.row(i).norm();
            radius += gi;
            g(i) = gi;
            if (i < dim - 1) sum += gi;
        }
        e = e * B;
        radius += e.norm();
        radius = 1.0 / radius;
        g = Bg * g;
        g = radius * g;
        for (i=0; i<dim; i++) temp_p.push_back(p0[i] + g(i));

        c = Point(dim, temp_p.begin(), temp_p.end());
        result.first = c;
        result.second = radius;

        return result;
    }


    /*
    // pick d+1 random vertices until they define a full dimensional simplex and then
    // compute the chebychev ball of that simplex
    std::pair<Point,NT> ComputeInnerBall() {

        std::vector<Point> verts(_d+1);
        std::vector<NT> vecp(_d);
        unsigned int vert_rand, pointer=0;
        unsigned int i,j;
        int m = num_of_vertices();
        std::vector<int> x_vec(_d);
        bool done=false;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_int_distribution<> uidist(1, m);

        std::pair<Point,NT> res;
        // while d+1 points do not define a full dimensional simplex repeat
        while(!done){
            pointer=0;
            x_vec.assign(_d+1,0);
            while(pointer!=_d+1) {
                vert_rand = uidist(rng);
                // Check if this vertex is selected first time
                if (std::find(x_vec.begin(), x_vec.begin() + pointer, vert_rand) == x_vec.begin() + pointer) {
                    x_vec[pointer] = vert_rand;
                    pointer++;
                }
            }

            for(i=0; i<(_d+1); i++){
                for (j=0; j<_d; j++) {
                    vecp[j] = V(x_vec[i] - 1, j);
                }
                verts[i] = Point(_d,vecp.begin(),vecp.end());
            }
            res=get_center_radius_inscribed_simplex(verts.begin(), verts.end(), done);
        }
        return res;
    }*/


    std::pair<Point,NT> ComputeInnerBall() {
        
        NT radius =  std::numeric_limits<NT>::max(), min_plus;
        Point center(_d);

        std::list<Point> randPoints;
        if (!get_points_for_rounding(randPoints)) {
            center = get_mean_of_vertices();
        } else {

            MT Ap(_d,randPoints.size());
            typename std::list<Point>::iterator rpit=randPoints.begin();

            unsigned int i, j = 0;
            for ( ; rpit!=randPoints.end(); rpit++, j++) {
                const NT* point_data = rpit->getCoefficients().data();

                for ( i=0; i < rpit->dimension(); i++){
                    Ap(i,j)=double(*point_data);
                    point_data++;
                }
            }
            MT Q(_d, _d);
            VT c2(_d);
            size_t w=1000;
            
            KhachiyanAlgo(Ap,0.01,w,Q,c2); // call Khachiyan algorithm

            //Get ellipsoid matrix and center as Eigen objects
            for(unsigned int i=0; i<_d; i++) center.set_coord(i, NT(c2(i)));
        }

        std::pair<NT,NT> res;
        Point v(_d);
        for (unsigned int i = 0; i < _d; ++i) {
            v.set_to_origin();
            v.set_coord(i, 1.0);
            res = intersect_double_line_Vpoly<NT>(V, center, v, row, colno);
            min_plus = std::min(res.first, -1.0*res.second);
            if (min_plus < radius) radius = min_plus;
        }

        radius = radius / std::sqrt(NT(_d));
        _inner_ball = std::pair<Point, NT> (center, radius);
        return std::pair<Point, NT> (center, radius);
    }


    // check if point p belongs to the convex hull of V-Polytope P
    int is_in(const Point &p, NT tol=NT(0)) const {
        if (memLP_Vpoly(V, p, conv_mem, colno_mem)){
            return -1;
        }
        return 0;
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(const Point &r, const Point &v) const {

        return intersect_double_line_Vpoly<NT>(V, r, v, row, colno);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(const Point &r, const Point &v, const VT &Ar,
            const VT &Av) const {
        return intersect_double_line_Vpoly<NT>(V, r, v,  row, colno);
    }

    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(const Point &r, const Point &v, const VT &Ar,
                                    const VT &Av, const NT &lambda_prev) const {

        return intersect_double_line_Vpoly<NT>(V, r, v,  row, colno);
    }


    std::pair<NT, int> line_positive_intersect(const Point &r, const Point &v) const {
        return std::pair<NT, int> (intersect_line_Vpoly(V, r, v, conv_comb, row, colno, false, false), 1);
    }

    std::pair<NT, int> line_positive_intersect(const Point &r, const Point &v, const VT &Ar,
                                               const VT &Av) const {
        return line_positive_intersect(r, v);
    }


    std::pair<NT, int> line_positive_intersect(const Point &r, const Point &v, const VT &Ar,
                                               const VT &Av, const NT &lambda_prev) const {
        return line_positive_intersect(r, v);
    }

    //-------------------------accelarated billiard--------------------------------//
    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(Point const& r,
                                                     Point const& v,
                                                     VT& Ar,
                                                     VT& Av,
                                                     update_parameters &params) const
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
                                               update_parameters &params) const
    {
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               update_parameters &params) const
    {
        return line_positive_intersect(r, v);
    }
    //------------------------------------------------------------------------------//


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(const Point &r,
                                          const unsigned int rand_coord,
                                          const VT &lamdas) const {
        Point v(_d);
        v.set_coord(rand_coord, 1.0);
        return intersect_double_line_Vpoly<NT>(V, r, v,  row, colno);
    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(const Point &r,
                                          const Point &r_prev,
                                          const unsigned int rand_coord,
                                          const unsigned int rand_coord_prev,
                                          const VT &lamdas) const {
        return line_intersect_coord(r, rand_coord, lamdas);
    }


    //------------------------------oracles for exponential sampling---------------//////

    // compute intersection points of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> quadratic_positive_intersect(Point const& r,
                                    Point const& v,
                                    VT const& Ac,
                                    NT const& T,
                                    VT& Ar,
                                    VT& Av,
                                    int& facet_prev) const
    {
        throw std::runtime_error("Quadratic polynomial trajectories are supported only for H-polytopes");
    }

    std::pair<NT, int> quadratic_positive_intersect(Point const& r,
                                    Point const& v,
                                    VT const& Ac,
                                    NT const& T,
                                    VT& Ar,
                                    VT& Av,
                                    NT const& lambda_prev,
                                    int& facet_prev) const
    {
        throw std::runtime_error("Quadratic polynomial trajectories are supported only for H-polytopes");
    }


    //------------oracle for exact hmc spherical gaussian sampling---------------//
    std::pair<NT, int> trigonometric_positive_intersect(Point const& r, Point const& v,
                                                      NT const& omega, int &facet_prev) const
    {
        return std::make_pair(0, 0);
    }


    // shift polytope by a point c
    void shift(const VT &c) {
        MT V2 = V.transpose().colwise() - c;
        V = V2.transpose();
    }


    // apply linear transformation, of square matrix T, to the V-Polytope
    void linear_transformIt(const MT &T) {
        MT V2 = T.inverse() * V.transpose();
        V = V2.transpose();
    }


    // consider an upper bound for the number of facets of a V-polytope
    // for each facet consider a lower bound for the distance from the origin
    // useful for CV algorithm to get the first gaussian
    std::vector<NT> get_dists(const NT &radius) const {
        std::vector <NT> res(upper_bound_of_hyperplanes(), radius);
        return res;
    }


    // in number_of_vertices<=20*dimension use the vertices for the rounding
    // otherwise you have to sample from the V-polytope
    template <class PointList>
    bool get_points_for_rounding (PointList &randPoints) {
        if (num_of_vertices()>20*_d) {
            return false;
        }
        unsigned int j;

        typename std::vector<NT>::iterator pointIt;
        for (int i=0; i<num_of_vertices(); i++) {
            Point p(V.row(i));
            randPoints.push_back(p);
        }
        return true;
    }

    void compute_reflection(Point &v, const Point &p, const int &facet) const {

        //compute_reflection(v, p, 0.0);

        int count = 0, outvert;
        MT Fmat2(_d,_d);
        for (int j = 0; j < num_of_vertices(); ++j) {
            if (*(conv_comb + j) > 0.0) {
                Fmat2.row(count) = V.row(j);
                count++;
            } else {
                outvert = j;
            }
        }

        VT a = Fmat2.colPivHouseholderQr().solve(VT::Ones(_d));
        if (a.dot(V.row(outvert)) > 1.0) a = -a;
        a /= a.norm();

        // compute reflection
        a *= (-2.0 * v.dot(a));
        v += a;
    }

    template <typename update_parameters>
    void compute_reflection(Point &v, const Point &p, update_parameters const& params) const {

        int count = 0, outvert;
        MT Fmat2(_d,_d);
        for (int j = 0; j < num_of_vertices(); ++j) {
            if (*(conv_comb + j) > 0.0) {
                Fmat2.row(count) = V.row(j);
                count++;
            } else {
                outvert = j;
            }
        }

        VT a = Fmat2.colPivHouseholderQr().solve(VT::Ones(_d));
        if (a.dot(V.row(outvert)) > 1.0) a *= -1.0;
        a /= a.norm();

        // compute reflection
        a *= (-2.0 * v.dot(a));
        v += a;
    }

    template <class bfunc, class NonLinearOracle>
    std::tuple<NT, Point, int> curve_intersect(
        NT t_prev,
        NT t0,
        NT eta,
        std::vector<Point> &coeffs,
        bfunc phi,
        bfunc grad_phi,
        NonLinearOracle &intersection_oracle,
        int ignore_facet=-1)
    {
        return intersection_oracle.apply(t_prev, t0, eta, V, *this, 
                                         coeffs, phi, grad_phi, ignore_facet);
    }


};

#endif
