// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HMC_HEUR_SAM_H
#define HMC_HEUR_SAM_H


template <class RNGType, class VT, class MT, typename NT>
void get_next_hmc_barrier(VT &x0, VT &v0, VT &s0, VT &sv0, MT &At, VT &b, MT &A, NT &a) {

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0, 1);
    boost::random::uniform_real_distribution<> urdist(0, 1);

    NT tpos = std::numeric_limits<NT>::max(), tminus = std::numeric_limits<NT>::lowest(), t1, t2;

    unsigned int m = s0.size();
    unsigned int d = x0.size();

    MT abcs(m, 3);
    abcs.col(2) = s0;
    abcs.col(1) = sv0;

    VT bt = b - abcs.col(2);
    for (int j = 0; j < m; ++j) bt(j) = a / bt(j);

    abcs.col(0) = (-0.5) * (At * bt);

    MT pinvAabcs(d,3);
    pinvAabcs.col(2) = x0;
    pinvAabcs.col(1) = v0;
    pinvAabcs.col(0) = abcs.col(0);

    abcs.col(0) = A * abcs.col(0);
    abcs.col(2) = abcs.col(2) - b;

    VT deltas = abcs.col(1).array()*abcs.col(1).array() - 4.0 * abcs.col(0).array() * abcs.col(2).array();
    for (int i = 0; i < m; ++i) {

        if (deltas(i) < 0.0) continue;

        t1 = (-abcs(i, 1) - std::sqrt(deltas(i))) / (2.0 * abcs(i, 0));
        if (t1 > 0 && t1 < tpos) tpos = t1;
        if (t1 < 0 && t1 > tminus) tminus = t1;

        t2 = (-abcs(i, 1) + std::sqrt(deltas(i))) / (2.0 * abcs(i, 0));
        if (t2 > 0 && t2 < tpos) tpos = t2;
        if (t2 < 0 && t2 > tminus) tminus = t2;

    }

    NT trand = urdist(rng) * (tpos - tminus) + tminus;
    VT tvec(3);
    tvec(0) = trand * trand;
    tvec(1) = trand;
    tvec(2) = 1.0;
    x0 = pinvAabcs * tvec;

}


template <class RNGType, class Polytope, class Point, class PointList, typename NT>
void hmc_logbarrier(Polytope &P, Point &p, PointList &randPoints, NT &a, int N) {

    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);

    MT A = P.get_mat();
    VT b = P.get_vec();
    unsigned int d = P.dimension();
    unsigned int m = A.rows();

    VT v0(d), x0(d);
    for (int i = 0; i < d; ++i) {
        v0(i) = rdist(rng);
        x0(i) = p[i];
    }

    VT s0 = A*x0;
    VT sv0 = A*v0;
    MT At = A.transpose();

    for (int i = 0; i < N; ++i) {

        get_next_hmc_barrier<RNGType>(x0, v0, s0, sv0, At, b, A, a);

        for (int i = 0; i < d; ++i) {
            v0(i) = rdist(rng);
            p.set_coord(i, x0(i) );
        }

        randPoints.push_back(p);
        s0 = A*x0;
        sv0 = A*v0;
    }

}



#endif
