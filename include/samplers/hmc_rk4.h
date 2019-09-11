// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HMC_RK4_H
#define HMC_RK4_H


template <class RNGType, class Polytope, class Point, class PointList, typename NT>
void hmc_logbarrier_rk2(Polytope &P, Point &p, int walk_step, PointList &randPoints, NT &a, int N,  NT radius = -1.0,
                        NT R = -1.0) {

    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::normal_distribution<> rdist(0, 1);

    MT A = P.get_mat(), At = A.transpose();
    unsigned int d = P.dimension(), m = A.rows();
    VT b = P.get_vec(), s0(m), Y1(d), Y2(d), Y051(d), m1i(d), m2i(d);
    NT sumh, h, T, r = 1.0, L = 1.0;
    const NT har = 0.1;

    if (R < 0.0) R = 1.0;
    if (radius > 0.0) {

        r = R * radius;
        A = A * (MT::Identity(d, d) * r);
        At = A.transpose();

    }

    for (int i = 0; i < d; ++i) Y1(i) = (1.0 / r) * p[i];

    for (int i = 0; i < N; ++i) {

        for (int l = 0; l < walk_step; ++l) {

            T = urdist(rng) * L;
            for (int i = 0; i < d; ++i) Y2(i) = rdist(rng);
            if (urdist(rng) > 0.5) Y2 = -Y2;

            sumh = 0.0;
            h = (T > har) ? har : har - T;

            while (sumh < T) {

                m1i = VT::Zero(d);
                //m2i = VT::Zero(d);

                for (int k = 0; k < 2; ++k) {
                    s0 = A * (Y1 + (0.5 * h) * m1i);
                    //(k < 1) ? A * Y1 :
                    m1i = (k < 1) ? Y2 : Y2 + (0.5 * h) * m2i;
                    for (int j = 0; j < m; ++j) s0(j) = a / (b(j) - s0(j));
                    m2i = -At * s0;
                }

                Y051 = Y1 + h * m1i;
                s0 = A * Y051;

                if ((s0 - b).maxCoeff() > 0.0) {
                    h = h * 0.25;
                    continue;
                }

                Y1 = Y051;
                Y2 = Y2 + h * m2i;
                sumh += h;
                h = (T - sumh > har) ? har : T - sumh;
            }
        }

        for (int k = 0; k < d; ++k) p.set_coord(k, r * Y1(k));
        randPoints.push_back(p);

    }
}


template <class RNGType, class Polytope, class Point, class PointList, typename NT>
void hmc_logbarrier_rk4(Polytope &P, Point &p, int walk_step, PointList &randPoints, NT &a, int N,  NT radius = -1.0,
                        NT R = -1.0) {

    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::normal_distribution<> rdist(0, 1);

    MT A = P.get_mat(), At = A.transpose();
    unsigned int d = P.dimension(), m = A.rows();
    VT b = P.get_vec(), s0(m), Y1(d), Y2(d), Y051(d), m1i(d), m2i(d), ms1(d), ms2(d);
    NT sumh, h, T, r = 1.0, L = 1.0;
    const NT har = 0.1;

    if (R < 0.0) R = 1.0;
    if (radius > 0.0) {

        r = R * radius;
        A = A * (MT::Identity(d, d) * r);
        At = A.transpose();

    }

    for (int i = 0; i < d; ++i) Y1(i) = (1.0 / r) * p[i];

    for (int i = 0; i < N; ++i) {

        for (int l = 0; l < walk_step; ++l) {

            T = urdist(rng) * L;
            for (int i = 0; i < d; ++i) Y2(i) = rdist(rng);
            if (urdist(rng) > 0.5) Y2 = -Y2;

            sumh = 0.0;
            h = (T > har) ? har : har - T;

            while (sumh < T) {

                ms1 = VT::Zero(d);
                ms2 = VT::Zero(d);
                m1i = VT::Zero(d);
                m2i = VT::Zero(d);

                for (int k = 0; k < 4; ++k) {
                    s0 = (k < 3) ? A * (Y1 + (0.5 * h) * m1i) : A * (Y1 + h * m1i);
                    m1i = (k < 3) ? Y2 + (0.5 * h) * m2i : Y2 + h * m2i;
                    for (int j = 0; j < m; ++j) s0(j) = a / (b(j) - s0(j));
                    m2i = -At * s0;
                    ms1 += (k > 0 && k < 3) ? 2.0 * m1i : m1i;
                    ms2 += (k > 0 && k < 3) ? 2.0 * m2i : m2i;
                }

                Y051 = Y1 + (h / 6.0) * ms1;
                s0 = A * Y051;

                if ((s0 - b).maxCoeff() > 0.0) {
                    h = h * 0.25;
                    continue;
                }

                Y1 = Y051;
                Y2 = Y2 + (h / 6.0) * ms2;
                sumh += h;
                h = (T - sumh > har) ? har : T - sumh;
            }
        }

        for (int k = 0; k < d; ++k) p.set_coord(k, r * Y1(k));
        randPoints.push_back(p);

    }
}


#endif
