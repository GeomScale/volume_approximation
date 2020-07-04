// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef GET_FULL_DIMENSIONAL_POLYTOPE
#define GET_FULL_DIMENSIONAL_POLYTOPE

/*template <typename MT>
MT get_Q(MT A){
    int n = A.cols();
    int m = A.rows();
    MT R_1(n,n);
    MT Q_1(m,n);
    for (int i = 0; i < n; ++i) {
        R_1(i,i) = A.col(i).norm();
        Q_1.col(i) = A.col(i) / R_1(i,i);
        for (int j = i+1; j < n; ++j) {
            R_1(i,j) = Q_1.col(i).transpose() * A.col(j);
            A.col(j) = A.col(j) - (Q_1.col(i)*R_1(i,j));
        }
    }
    return Q_1;
}*/



template <typename H_polytope, typename MT, typename VT>
std::pair<H_polytope, std::pair<MT, VT> > get_full_dimensional_polytope(MT A, VT b, MT Aeq, VT beq)
{
    typedef typename H_polytope::NT NT;

    VT p = Aeq.colPivHouseholderQr().solve(beq);

    //std::cout<<"p = "<<p.transpose()<<std::endl;

    //Eigen::FullPivLU<MT> lu(Aeq);
    //MT N = lu.kernel();

    Eigen::CompleteOrthogonalDecomposition<MT> cod;
    cod.compute(Aeq);

    Eigen::ColPivHouseholderQR<MT> qrdecomp(Aeq.transpose());
    MT N = qrdecomp.householderQ();
    //std::cout << "Q:\n" << Q << std::endl;
 
    // Find URV^T
    //MT V = cod.matrixZ().transpose();
    //MT N = V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());
    //MT P = cod.colsPermutation();
    //N = P * N; // Unpermute the columns

    b = b - A * p;
    A = A * N;

    H_polytope HP;
    HP.init(A.cols(), A, b);

    return std::pair<H_polytope, std::pair<MT, VT> >(HP, std::pair<MT,VT>(N, p));

}

#endif
