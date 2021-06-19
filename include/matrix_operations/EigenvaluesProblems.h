// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_EIGENVALUESPROBLEMS_H
#define VOLESTI_EIGENVALUESPROBLEMS_H

/// Uncomment the solver the function minPosGeneralizedEigenvalue uses
/// Eigen solver for generalized eigenvalue problem
//#define EIGEN_EIGENVALUES_SOLVER
/// Spectra standard eigenvalue problem
#define SPECTRA_EIGENVALUES_SOLVER
/// ARPACK++ standard eigenvalues solver
//#define ARPACK_EIGENVALUES_SOLVER

//#include <../../external/arpack++/include/arssym.h>
#include <../../external/Spectra/include/Spectra/SymEigsSolver.h>
#include "DenseProductMatrix.h"
#include "EigenDenseMatrix.h"

#include "../../external/Spectra/include/Spectra/SymGEigsSolver.h"
#include "../../external/Spectra/include/Spectra/GenEigsSolver.h"
//#include "../../external/arpack++/include/arsnsym.h"

/// Solve eigenvalues problems
/// \tparam NT Numeric Type
/// \tparam MT Matrix Type
/// \tparam VT Vector Type
template<typename NT, typename MT, typename VT>
class EigenvaluesProblems {

};


/// A specialization of the template class EigenvaluesProblems for dense Eigen matrices and vectors.
/// \tparam NT
template<typename NT>
class EigenvaluesProblems<NT, Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>, Eigen::Matrix<NT,Eigen::Dynamic,1> > {
public:
    /// The type for Eigen Matrix
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    /// The type for Eigen vector
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    /// The type of a complex Eigen Vector for handling eigenvectors
#if defined(EIGEN_EIGENVALUES_SOLVER) || defined (SPECTRA_EIGENVALUES_SOLVER)
    typedef typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType CVT;
#elif defined(ARPACK_EIGENVALUES_SOLVER)
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> CVT;
#endif

    /// The type of a pair of NT
    typedef std::pair<NT, NT> NTpair;


    /// Find the smallest eigenvalue of M
    /// \param M a symmetric matrix
    /// \return smallest eigenvalue
    NT findSymEigenvalue(MT const & M) {
        EigenDenseMatrix<NT> _M(&M);

//#define NOT_WORKING
#ifdef NOT_WORKING
        // Creating an eigenvalue problem and defining what we need:
        // the smallest eigenvalue of M.
        ARNonSymStdEig<NT, EigenDenseMatrix<NT> >
                dprob(M.cols(), 1, &_M, &EigenDenseMatrix<NT>::MultMv, std::string ("LR"), 8, 0.0, 100*15);

        // compute
        if (dprob.FindEigenvectors() == 0) {
            std::cout << "Failed in findSymEigenvalue\n";
            // if failed with default (and fast) parameters, try with stable (and slow)
            dprob.ChangeNcv(M.cols()/10);
            if (dprob.FindEigenvectors() == 0) {
                std::cout << "\tFailed Again\n";
                return NT(0);
            }
        }

        if (!dprob.EigenvaluesFound()) {
            // if failed to find eigenvalues
            return NT(0);
        }

        // retrieve eigenvalue of the original system
        return dprob.EigenvalueReal(0);
#elif defined(SPECTRA)
        // This parameter is for Spectra. It must be larger than #(requested eigenvalues) + 2
        // and smaller than the size of matrix;
        int ncv = M.cols()/10 + 5;
        if (ncv > M.cols()) ncv = M.cols();

        Spectra::SymEigsSolver<NT, Spectra::LARGEST_ALGE, EigenDenseMatrix<NT> > eigs(&_M, 1, ncv);
        // compute
        eigs.init();
        eigs.compute(50000);
        if(eigs.info() == Spectra::SUCCESSFUL) {
            return eigs.eigenvalues()(0);
        }
        else {
            std::cout << "Spectra failed\n";
            return NT(0);
        }
#else
        Eigen::SelfAdjointEigenSolver<MT> solver;
        solver.compute(M, Eigen::EigenvaluesOnly);
//        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();
//        NT max = eivals(0).real();
//
//        for (int i = 1; i < eivals.rows(); i++)
//            if (eivals(i).real() > max)
//                max = eivals(i).real();

        return solver.eigenvalues().maxCoeff();
#endif
    }

    /// Find the minimum positive and maximum negative eigenvalues of the generalized eigenvalue
    /// problem A + lB, where A, B symmetric and A negative definite.
    /// \param[in] A Input matrix
    /// \param[in] B Input matrix
    /// \return The pair (minimum positive, maximum negative) of eigenvalues
    NTpair symGeneralizedProblem(MT const & A, MT const & B) {

        int matrixDim = A.rows();

        // Spectra solves Xv=lYv, where Y positive definite
        // Set X = B, Y=-A. Then, the eigenvalues we want are the minimum negative
        // and maximum positive eigenvalues of Xv=lYv.

        // Construct matrix operation object using the wrapper classes provided by Spectra
        Spectra::DenseSymMatProd<NT> op(B);
        Spectra::DenseCholesky<NT> Bop(-A);

        // Construct generalized eigen solver object
        // requesting the minmum negative and largest positive eigenvalues
        Spectra::SymGEigsSolver<NT, Spectra::BOTH_ENDS, Spectra::DenseSymMatProd<NT>, Spectra::DenseCholesky<NT>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 2, 5 < matrixDim ? 5 : matrixDim);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();

        // Retrieve results
        if (geigs.info() != Spectra::SUCCESSFUL)
            return {NT(0), NT(0)};

        Eigen::VectorXd evalues;
        double lambdaMinPositive, lambdaMaxNegative;

        evalues = geigs.eigenvalues();

        // get the eigenvalues of the original problem
        lambdaMinPositive = 1 / evalues(0);
        lambdaMaxNegative = 1 / evalues(1);

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    /// Finds the minimum positive real eigenvalue of the generalized eigenvalue problem A + lB and
    /// the corresponding eigenvector.
    /// If the macro EIGEN_EIGENVALUES_SOLVER is defined, the Generalized Solver of Eigen is used.
    /// Otherwise, we transform the generalized to a standard eigenvalue problem and use Spectra.
    /// Warning: With Spectra we might get a value smaller than the minimum positive real eigenvalue (the real part
    /// of a complex eigenvalue).
    /// No restriction on the matrices!
    /// \param[in] A Input matrix
    /// \param[in] B Input matrix
    /// \param[out] eigenvector The eigenvector corresponding to the minimum positive eigenvalue
    /// \return The minimum positive eigenvalue
    NT minPosGeneralizedEigenvalue(MT const & A, MT const & B, CVT& eigenvector) {
        NT lambdaMinPositive = std::numeric_limits<NT>::max();

#if defined(EIGEN_EIGENVALUES_SOLVER)
        // use the Generalized eigenvalue solver of Eigen

        // compute generalized eigenvalues with Eigen solver
        Eigen::GeneralizedEigenSolver<MT> ges(A, -B);

        // retrieve minimum positive eigenvalue
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();
        int index = 0;

        for (int i = 0; i < alphas.rows(); i++) {

            if (betas(i) == 0 || alphas(i).imag() != 0)
                continue;

            double lambda = alphas(i).real() / betas(i);
            if (lambda > 0 && lambda < lambdaMinPositive) {
                lambdaMinPositive = lambda;
                index = i;
            }
        }

        // retrieve corresponding eigenvector
        eigenvector = ges.eigenvectors().col(index);
#elif defined(SPECTRA_EIGENVALUES_SOLVER)
        // Transform the problem to a standard eigenvalue problem and use the general eigenvalue solver of Spectra

        // This makes the transformation to standard eigenvalue problem. See class for more info.
        // We have the generalized problem  A + lB, or Av = -lBv
        // This class computes the matrix product vector Mv, where M = -B * A^[-1]
        MT _B = -1 * B; // TODO avoid this allocation
        DenseProductMatrix<NT> M(&_B, &A);

        // This parameter is for Spectra. It must be larger than #(requested eigenvalues) + 2
        // and smaller than the size of matrix;
        int ncv = 3;

        // Prepare to solve Mx = (1/l)x
        // we want the smallest positive eigenvalue in the original problem,
        // so in this the largest positive eigenvalue;
        Spectra::GenEigsSolver<NT, Spectra::LARGEST_REAL, DenseProductMatrix<NT> > eigs(&M, 1, ncv);

        // compute
        eigs.init();
        eigs.compute();

        //retrieve result and invert to get required eigenvalue of the original problem
        if (eigs.info() != Spectra::SUCCESSFUL) {
            eigenvector.setZero(A.rows());
            return NT(0);
        }

        lambdaMinPositive = 1/((eigs.eigenvalues())(0).real());

        // retrieve corresponding eigenvector
        int matrixDim = A.rows();
        eigenvector.resize(matrixDim);
        for (int i = 0; i < matrixDim; i++)
            eigenvector(i) =  (eigs.eigenvectors()).col(0)(i);

#elif defined(ARPACK_EIGENVALUES_SOLVER)
        // Transform the problem to a standard eigenvalue problem and use the general eigenvalue solver of ARPACK++

        // This makes the transformation to standard eigenvalue problem. See class for more info.
        // We have the generalized problem  A + lB, or Av = -lBv
        // This class computes the matrix product vector Mv, where M = -B * A^[-1]
        MT _B = -1 * B; // TODO avoid this allocation
        DenseProductMatrix<NT> M(&_B, &A);

        // Creating an eigenvalue problem and defining what we need:
        // the  eigenvector of A with largest real.
        ARNonSymStdEig<NT, DenseProductMatrix<NT> >

        dprob(A.cols(), 1, &M, &DenseProductMatrix<NT>::MultMv, std::string ("LR"), 8<A.rows() ? 8 : A.rows(), 0.000);//, 100*3);

        // compute
        if (dprob.FindEigenvectors() == 0) {
            std::cout << "Failed\n";
            // if failed with default (and fast) parameters, try with stable (and slow)
            dprob.ChangeNcv(A.cols()/10);
            if (dprob.FindEigenvectors() == 0) {
                std::cout << "\tFailed Again\n";
                return NT(0);
            }
        }


        // allocate memory for the eigenvector here
        eigenvector.setZero(A.rows());

        if (!dprob.EigenvaluesFound()) {
            // if failed to find eigenvalues
            return NT(0);
        }

        // retrieve eigenvalue of the original system
        lambdaMinPositive = 1/dprob.EigenvalueReal(0);

        eigenvector.setZero(A.rows());
        if (dprob.EigenvectorsFound()) {
            //retrieve corresponding eigenvector
            for (int i=0 ;i<A.rows() ; i++)
                eigenvector(i) = dprob.EigenvectorReal(0, i);
        }


#endif
//        std::cout << lambdaMinPositive << " " << eigenvector.transpose() << "\n";fflush(stdout);
        return lambdaMinPositive;
    }

    /// Find the minimum positive and maximum negative eigenvalues of the generalized eigenvalue
    /// problem A + lB, where A, B symmetric and A negative definite.
    /// \param[in] A Input matrix
    /// \param[in] B Input matrix
    /// \return The pair (minimum positive, maximum negative) of eigenvalues
    NT minPosLinearEigenvalue(MT const & A, MT const & B, VT &eigvec) 
    {
        int matrixDim = A.rows();
        double lambdaMinPositive;

        Spectra::DenseSymMatProd<NT> op(B);
        Spectra::DenseCholesky<NT> Bop(-A);
        //std::cout<<"A = "<<Eigen::MatrixXd(A)<<"\n\n"<<std::endl;
        //std::cout<<"B = "<<Eigen::MatrixXd(B)<<"\n\n"<<std::endl;

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
        Spectra::SymGEigsSolver<NT, Spectra::LARGEST_ALGE,  Spectra::DenseSymMatProd<NT>, Spectra::DenseCholesky<NT>, Spectra::GEIGS_CHOLESKY> 
            geigs(&op, &Bop, 1, 15 < matrixDim ? 15 : matrixDim);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();

        // Retrieve results
        VT evalues;
        //Eigen::MatrixXd evecs;

        if (geigs.info() == Spectra::SUCCESSFUL) {
            evalues = geigs.eigenvalues();
            eigvec = geigs.eigenvectors().col(0);
        }

        lambdaMinPositive = 1 / evalues(0);

        //std::cout<<"lambdaMinPositive = "<<lambdaMinPositive<<std::endl;
        //std::cout<<"eigvec = "<<eigvec.transpose()<<"\n\n"<<std::endl;

        //std::cout<<"lBx = "<<lambdaMinPositive*(B.template selfadjointView< Eigen::Lower >()*eigvec).transpose()<<"\n"<<std::endl;
        //std::cout<<"Ax = "<<(((-A).template selfadjointView< Eigen::Lower >()*eigvec).transpose())<<"\n\n"<<std::endl;

        return lambdaMinPositive;
    }

    /// Transform the quadratic eigenvalue problem \[At^2 + Bt + c\] to
    /// the generalized eigenvalue problem X+lY.
    /// If the updateOnly flag is false, compute matrices X,Y from scratch;
    /// otherwise update them.
    /// \param[in] A
    /// \param[in] B
    /// \param[in] C
    /// \param[in, out] X
    /// \param[in, out] Y
    /// \param[in, out] updateOnly True if X,Y were previously computed and only B,C changed
    void linearization(const MT &A, const MT &B, const MT &C, MT &X, MT &Y, bool &updateOnly) {
        unsigned int matrixDim = A.rows();

        // check if the matrices X,Y are computed.
        //if yes, update them; otherwise compute them from scratch
        if (!updateOnly) {
            X.resize(2 * matrixDim, 2 * matrixDim);
            Y.resize(2 * matrixDim, 2 * matrixDim);

            Y.block(matrixDim, matrixDim, matrixDim, matrixDim) = -1 * C;
            Y.block(0, matrixDim, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
            Y.block(matrixDim, 0, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
            Y.block(0, 0, matrixDim, matrixDim) = A;

            X.block(0, matrixDim, matrixDim, matrixDim) = C;
            X.block(0, 0, matrixDim, matrixDim) = B;
            X.block(matrixDim, 0, matrixDim, matrixDim) = C;
            X.block(matrixDim, matrixDim, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
        } else {
            Y.block(matrixDim, matrixDim, matrixDim, matrixDim) = -1 * C;

            X.block(0, matrixDim, matrixDim, matrixDim) = C;
            X.block(0, 0, matrixDim, matrixDim) = B;
            X.block(matrixDim, 0, matrixDim, matrixDim) = C;
        }
    }

    /// Find the minimum positive real eigenvalue of the quadratic eigenvalue problem \[At^2 + Bt + c\].
    /// First transform it to the generalized eigenvalue problem X+lY.
    /// If the updateOnly flag is false, compute matrices X,Y from scratch;
    /// otherwise only update them.
    /// \param[in] A Input matrix
    /// \param[in] B Input matrix
    /// \param[in] C Input matrix
    /// \param[in, out] X
    /// \param[in, out] Y
    /// \param[out] eigenvector The eigenvector corresponding to the minimum positive eigenvalue
    /// \param[in, out] updateOnly True if X,Y were previously computed and only B,C changed
    /// \return Minimum positive eigenvalue
    NT
    minPosQuadraticEigenvalue(MT const & A, MT const &B, MT const &C, MT &X, MT &Y, VT &eigenvector, bool &updateOnly) {
        // perform linearization and create generalized eigenvalue problem X+lY
        linearization(A, B, C, X, Y, updateOnly);

        // solve generalized problem
        CVT eivector;
        NT lambdaMinPositive = minPosGeneralizedEigenvalue(X, Y, eivector);

        if (lambdaMinPositive == 0)
            return 0;

        int matrixDim = A.rows();

        // the eivector has dimension 2*matrixDim
        // while the eigenvector of the original problem has dimension matrixDim
        // retrieve the eigenvector by keeping only #matrixDim coordinates.
        eigenvector.resize(matrixDim);

#if defined(EIGEN_EIGENVALUES_SOLVER) || defined (SPECTRA_EIGENVALUES_SOLVER)
        for (int i = 0; i < matrixDim; i++)
            eigenvector(i) =  eivector(matrixDim + i).real();
#elif defined(ARPACK_EIGENVALUES_SOLVER)
        for (int i = 0; i < matrixDim; i++)
            eigenvector(i) =  eivector(matrixDim + i);
#endif

        return lambdaMinPositive;
    }
};

#endif //VOLESTI_EIGENVALUESPROBLEMS_H
