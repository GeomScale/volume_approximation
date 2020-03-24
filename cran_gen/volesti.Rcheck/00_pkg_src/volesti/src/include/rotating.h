// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2019 Apostolos Chalkis


// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ROTATING_H
#define ROTATING_H



template <typename MT, typename Polytope>
MT rotating(Polytope &P){

    typedef boost::mt19937    RNGType;
    //typedef typename Polytope::MT 	MT;

    boost::random::uniform_real_distribution<> urdist(-1.0, 1.0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    unsigned int n = P.dimension();

    // pick a random rotation
    MT R(n,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            R(i,j) = urdist(rng);
        }
    }

    Eigen::JacobiSVD<MT> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // apply rotation to the polytope P
    P.linear_transformIt(svd.matrixU());

    return svd.matrixU().inverse();
}

#endif
