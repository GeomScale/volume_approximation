// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SPECTRAHEDRON_NEW_H
#define VOLESTI_SPECTRAHEDRON_NEW_H

#include <boost/format.hpp>

#include "../../cartesian_geom/cartesian_kernel.h"
#include "newLMI.h"
#include "chrono"
#include "generators/boost_random_number_generator.hpp"
#include "../../sampling/sphere.hpp"


template <typename MT>
struct PrecomputedValues {
    
};

/// Among successive calls of this class methods, we may need to pass data
/// from one call to the next, to avoid repeating computations, or to efficiently update values
/// Warning: this struct assists in many methods; perhaps for different methods use different instances
template <typename NT>
struct PrecomputedValues<Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>> {

    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    /// These flags indicate whether the corresponding matrices are computed
    /// if yes, we can use them and not compute them fro scratch
    bool computed_A = false;
    bool computed_C = false;
    bool computed_XY = false;

    /// The matrices the method positiveIntersection receives from its previous call
    /// if the flag first_positive_intersection is true.
    /// Matrix A is also used in coordinateIntersection
    MT A, B, C;//, X, Y;

    /// In method positive_intersect, the distance we are computing corresponds
    /// to the minimum positive eigenvalue of a quadratic eigenvalue problem.
    /// This will hold the eigenvector for that eigenvalue
    VT eigenvector;

    /// Sets all flags to false
    void resetFlags() {
        computed_XY = computed_C = computed_A = false;
    }

    void set_mat_size(int const& m) 
    {
        //matrix_initializer<MT>::initialize(A, B, C, X, Y, m);

        A.setZero(m, m);
        B.setZero(m, m);
        C.setZero(m, m);

        eigenvector.setZero(m);

        //X.setZero(2*m, 2*m);
        //Y.setZero(2*m, 2*m);
    }
};






/// This class manipulates a spectrahedron, described by a LMI
/// \tparam NT Numeric Type
/// \tparam MT Matrix Type
/// \tparam VT Vector Type
template<typename NTT, typename MTT, typename VTT>
class Spectrahedron {
public:

    /// The numeric/matrix/vector types we use
    typedef NTT NT;
    typedef MTT MT;
    typedef VTT VT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point PointType;
    typedef PrecomputedValues<MT> PrecomputedValues2;

    /// The type of a pair of NT
    typedef std::pair<NT, NT> pairNT;
    PrecomputedValues2 precomputedValues;
    EigenvaluesProblems<NT, MT, VT> EigenvaluesProblem;

    /// The dimension of the spectrahedron
    unsigned int d;
    VT grad;

    /// The linear matrix inequality that describes the spectrahedron
    LMI<NT, MT, VT> lmi;

    Spectrahedron() {}

    /// Creates a spectrahedron
    /// \param[in] lmi The linear matrix inequality that describes the spectrahedron
    Spectrahedron(const LMI<NT, MT, VT>& lmi) : lmi(lmi) {
        d = lmi.dimension();
        grad.setZero(d);
        precomputedValues.resetFlags();
        precomputedValues.set_mat_size(lmi.sizeOfMatrices());
    }

    MT get_mat() const
    {
        return MT::Identity(lmi.dimension(), lmi.dimension());
    }

    int num_of_hyperplanes() const
    {
        return 0;
    }

    void resetFlags() {
        precomputedValues.resetFlags();
    }


    /// Construct the quadratic eigenvalue problem \[At^2 + Bt + C \] for positive_intersect.
    /// A = lmi(c) - A0, B = lmi(b) - A0 and C = lmi(c).
    /// \param[in] a Input vector
    /// \param[in] b Input vector
    /// \param[in] c Input vector
    /// \param[in, out] precomputedValues Holds matrices A, C
    void createMatricesForPositiveIntersection(const VT& p, const VT& v) {
        // check if matrices A, C are ready
        // if not compute them
        if (!precomputedValues.computed_C) {
            //std::cout<<"to evaluate lmi"<<std::endl;
            lmi.evaluate(p, precomputedValues.C);
            //std::cout<<"lmi evaluated"<<std::endl;
        }

        //if (!precomputedValues.computed_C) {
        //    lmi.evaluate(c, precomputedValues.C);
        //}

        // compute Matrix B
        lmi.evaluateWithoutA0(v, precomputedValues.B);
    }

    NT positiveIntersection(VT const & p, VT const & v) 
    {
        //unsigned int matrixDim = lmi.sizeOfMatrices();

        // create matrices A, B, C
        //std::cout<<"to create matrices"<<std::endl;
        createMatricesForPositiveIntersection(p, v);
        //std::cout<<"matrices created"<<std::endl;

        // get the minimum positive eigenvalue of At^2 + Bt + C
        
        //#if defined(SPARSE_PROBLEM)
            //SparseEigenvaluesProblems<NT, MT, VT> EigenvaluesProblem;
        //#elif defined(DENSE_PROBLEM)
        //    EigenvaluesProblems<NT, MT, VT> EigenvaluesProblem;
        //#endif
        //std::cout<<"computing eigenvalue..."<<std::endl;
        NT distance = EigenvaluesProblem.minPosLinearEigenvalue(precomputedValues.C, precomputedValues.B,
                                                                            precomputedValues.eigenvector);
        //std::cout<<"distance = "<<distance<<std::endl;

        return distance;
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType& r,
                                                     PointType& v,
                                                     VT&,
                                                     VT& ,
                                                     NT&,
                                                     update_parameters&) 
    {
        NT pos_inter = positiveIntersection(r.getCoefficients(), v.getCoefficients());
        return std::pair<NT, int> (pos_inter, -1);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType& r,
                                               PointType& v,
                                               VT& ,
                                               VT& ,
                                               NT &,
                                               MT& ,
                                               update_parameters& ) 
    {
        NT pos_inter = positiveIntersection(r.getCoefficients(), v.getCoefficients());
        return std::pair<NT, int> (pos_inter, -1);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(PointType& r,
                                                     PointType& v,
                                                     VT& ,
                                                     VT& ,
                                                     update_parameters&)
    {
        NT pos_inter = positiveIntersection(r.getCoefficients(), v.getCoefficients());
        return std::pair<NT, int> (pos_inter, -1);
    }

    void update_position_internal(NT &t){
        precomputedValues.C += t * precomputedValues.B;
        precomputedValues.computed_C = true;
    }

    /// Computes the reflected direction at a point on the boundary of the spectrahedron.
    /// \param[in] point A point on the boundary of the spectrahedron
    /// \param[in] incomingDirection The direction of the trajectory as it hits the boundary
    /// \param[out] reflectedDirection The reflected direction
    /// \param[in] precomputedValues Must contain an eigenvalue needed to compute the reflection
    void computeReflection(VT const & point, VT & incomingDirection, VT& reflectedDirection) {

        // get the gradient of the determinant of the lmi at point
        //VT grad;
        lmi.normalizedDeterminantGradient(point, precomputedValues.eigenvector, grad);

        // compute reflected direction
        // if v is original direction and s the surface normal,
        // reflected direction = v - 2 <v,s>*s

        NT dot = 2 * incomingDirection.dot(grad);
        reflectedDirection = incomingDirection - dot * grad;
    }

    template <typename update_parameters>
    void compute_reflection(PointType &v, const PointType &r, update_parameters& ) {

        lmi.normalizedDeterminantGradient(r.getCoefficients(), precomputedValues.eigenvector, grad);

        // compute reflected direction
        // if v is original direction and s the surface normal,
        // reflected direction = v - 2 <v,s>*s

        NT dot = 2 * v.dot(grad);
        v += -dot * PointType(grad);
    }

    /// \return The dimension of the spectrahedron
    unsigned int dimension() const {
        return d;
    }

    /// \return The LMI describing this spectrahedron
    LMI<NT, MT, VT> getLMI() const {
        return lmi;
    }



    /// Find out is lmi(current position) = mat is in the exterior of the spectrahedron
    /// \param mat a matrix where mat = lmi(current position)
    /// \return true if position is outside the spectrahedron
    bool isExterior(MT const & mat) {
        return !lmi.isNegativeSemidefinite(mat);
    }

    /// Find out is pos is in the exterior of the spectrahedron
    /// \param pos a vector
    /// \return true if pos is outside the spectrahedron
    bool isExterior(VT &pos) {
        return !lmi.isNegativeSemidefinite(pos);
    }
};

#endif //VOLESTI_SPECTRAHEDRON_H
