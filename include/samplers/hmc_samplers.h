// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HMC_SAMPLERS_H
#define HMC_SAMPLERS_H

#include <cmath>

template <typename T, typename T2>
T extract(const T2& full, const T& ind)
{
    int num_indices = ind.innerSize();
    T target(num_indices);
    for (int i = 0; i < num_indices; i++)
    {
        target[i] = full[ind[i]];
    }
    return target;
}

template <class RNGType, class Polytope, class Point, class PointList, typename NT>
void hmc_logbarrier(Polytope &P, Point &p, PointList &randPoints, NT &a, int n, int N) {

    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);

    MT A = P.get_mat();
    VT b = P.get_vec();
    unsigned int d = P.dimension();
    unsigned int m = A.nrows();

    VT v0(d), x0(d);
    for (int i = 0; i < d; ++i) {
        v0(i) = rdist(rng);
        x0(i) = p[i];
    }

    s0 = A*x0;
    sv0 = A*v0;
    MT M = - A*A.transpose();

    VT c(3);
    cj(0) = 0.0; cj(1) = 0.587785; cj(2) = 0.951056; //Chebyshev nodes

    MT AA = MT::Ones(n+1, n+1);
    MT T = MT::Zero(n+1,n+1);
    VT S = VT::Zero(n+1);

    for (int j = 1; j < n+1; ++j) {
        AA.col(j) = cj.pow(j);
        S(j) = NT(j)*cj,pow(j-1);
        if (j>1) T.col(j) = NT(j)*NT(j-1)*cj.pow(j-2);
    }
    AAinv = A.inverse();
    T = T*AAinv;
    S = S*AAinv;
    MT pinvA = A.completeOrthogonalDecomposition().pseudoInverse();

    for (int i = 0; i < N; ++i) {

        get_next_hmc_logbarrier<RNGType>(x0, T, S, M, s0, sv0, cj, pinvA, A, b, a);

        for (int i = 0; i < d; ++i) {
            v0(i) = rdist(rng);
            p[i] = x0(i);
        }
        randPoints.push_back(p);
        s0 = A*x0;
        sv0 = A*v0;
    }

}


template <class RNGType, class VT, class MT, typename NT>
void get_next_hmc_logbarrier(VT &x0, MT &T, VT &S, MT &M, VT &s0, VT &sv0, VT &cj, MT pinvA, MT &A, VT &b, NT &a) {

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<> urdist(0,1);

    unsigned int m = s0.size();
    unsigned int d = x0.size();
    unsigned int n = cj.size()-1;

    NT tpos = std::numeric_limits<NT>::max(), tminus = std::numeric_limits<NT>::lowest(), g, Delta, t1, t2;

    MT J = MT::Zero(m*(n+1), m*(n+1));
    VT F = VT::Zero(M*(n+1));
    VT s1(m*(n+1));
    for (int j = 0; j < m*(n+1); ++j) s1 = rdist(rng);

    VT s2 = 10*s1;

    for (int i = 0; i < m; ++i) {
        J(m * (n - 1) + i, i * (n + 1)) = 1.0;
        for (int j = 0; j < n + 1; ++j) J(m * n + i, (i - 1) * (n + 1) + j) = S(1, j);


        for (int j = 1; j < n - 1; ++j) {
            for (int k = 0; k < n + 1; ++k) J((i - 1) * (n - 1) + j - 1, (i - 1) * (n + 1) + k) = T(j, k);

            for (int k = 0; k < m; ++k) vec(k) = (k - 1) * (n + 1) + j;
        }
    }

    NT pre_dif = (s1-s2).abs().MaxCoeff();
    int count = 0, iter =0;
    VT Jirows(n+1);//=(1:(n+1))-1;

    for (int i = 0; i < n+1; ++i) {
        Jirows(i) = i-1;
    }

    while ((s1-s2).abs().MaxCoeff()>=0.000001) {

        iter++;
        for (int i = 0; i < m; ++i) {

            F(m * (n - 1) + i) = s1((i - 1) * (n + 1)) - s0(i);
            F(m * n + i) = S.dot(s1.segment(((i - 1) * (n + 1)), ((i - 1) * (n + 1) + n))) - sv0(i);

            for (int j = 1; j < n - 1; ++j) {
                J((i - 1) * (n - 1) + j - 1, (i - 1) * (n + 1) + j) = T(j, j) -
                        M(i, i) * a * (1 / ((b(i) - s1((i - 1) * (n + 1) + j)) * (b(i) - s1((i - 1) * (n + 1) + j))));
                F((i - 1) * (n - 1) + j - 1) =
                        T.row(j).dot(s1.segment(((i - 1) * (n + 1)), ((i - 1) * (n + 1) + n))) -
                        a * M.row(i).dot(b - extract(s1, vec).pow(-1.0));
            }
        }

        update_svals(s1, s2, J, F, Jirows, n, m);

        if (std::abs(pre_dif-(s1-s2).abs().MaxCoeff())<0.000001 && (s1-s2).abs().MaxCoeff()>0.0005) {

            for (int j = 0; j < m*(n+1); ++j) s1 = rdist(rng);

            VT s2 = 10*s1;
            count++;
            iter = 0;
        }
        pre_dif = (s1-s2).abs().MaxCoeff();

    }

    MT S(m, n+1);
    for (int i = 0; i < m; ++i) S.row(i) = s2.segment((i-1)*(n+1) , (i-1)*(n+1) + n);

    MT X = pinvA * S;
    MT abcs(d, n+1);
    MT aa;
    aa << cj(1)*cj(1), cj(1),
          cj(2)*cj(2), cj(2);
    MT bb(2);

    for (int i = 0; i < d; ++i) {
        g = X(i,0);
        bb(0)=X(i,1) - g; b(1) = X(i,2) - g;
        abcs.row(i) = aa.colPivHouseholderQr().solve(bb);
    }
    MT abc = A*abcs;

    for (int i = 0; i < m; ++i) {

        abc(i,3) = abc(i,2) - b(i);
        Delta = abc(i,1)*abc(i,1)-4.0*abc(i,0)*abc(i,2);
        if (Delta < 0.0) continue;

        t1 = (-abc(i,1)-std::sqrt(Delta))/(2.0*abc(i,0));
        if (t1>0.0 && t1<tpos) tpos = t1;
        if (t1<0.0 && t1>tminus) tminus = t1;

        t2 = (-abc(i,1)+sqrt(Delta))/(2.0*abc(i,0));
        if (t2>0.0 && t2<tpos) tpos = t2;
        if (t2<0.0 && t2>tminus) tminus = t2;

    }

    NT trand = urdist(rng)*(tpos - tminus) + tminus;
    VT tvec(3);
    tvec(0) = trand*trand; tvec(1) = trand; tvec(2) = 1.0;
    x0 = abcs * tvec;
    
}


// Newton method
template <class VT, class MT>
void update_svals(VT &s1, VT &s2, MT &J, VT &F, VT &Jirows, int n, int m) {


    MT Ji(n+1,n+1);
    VT Fi(n+1);

    VT stemp = s2;

    for (int i = 0; i < m; ++i) {

        for (int j = 0; j < n+1; ++j) {

            Ji.row(j) = J.row(Jirows(j) * m + i).segment(((i - 1) * (n + 1)), ((i - 1) * (n + 1) + n));
            Fi(j) = F(Jirows(j) * m + i);

        }

        s2.segment(((i-1)*(n+1)) , ((i-1)*(n+1)+n)) = s1.segment(((i-1)*(n+1)) , ((i-1)*(n+1)+n))
                - Ji.colPivHouseholderQr().solve(Fi);

    }
    s1 = stemp;

}


#endif
