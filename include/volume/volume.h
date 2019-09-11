// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLUME_H
#define VOLUME_H

#include <iterator>
//#include <fstream>
#include <vector>
#include <list>
//#include <algorithm>
#include <math.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "vars.h"
#include "polytopes.h"
//#include "ellipsoids.h"
#include "ballintersectconvex.h"
#include "vpolyintersectvpoly.h"
#include "samplers.h"
#include "rounding.h"
#include "rotating.h"
#include "gaussian_samplers.h"
#include "gaussian_annealing.h"
//#include "sample_only.h"
#include "misc.h"
#include "linear_extensions.h"
//#include "polytope_generators.h"
//#include "exact_vols.h"
//#include "simplex_samplers.h"
//#include "copulas.h"


template <class Polytope, class Parameters, class Point, typename NT>
NT volume(Polytope &P,
                  Parameters &var,  // constans for volume
                  std::pair<Point,NT> InnerBall)  //Chebychev ball
{
 
    typedef Ball<Point> Ball;
    typedef BallIntersectPolytope<Polytope,Ball> BallPoly;
    typedef typename Parameters::RNGType RNGType;
    typedef typename Polytope::VT VT;

    bool round = var.round;
    bool print = var.verbose;
    bool rand_only = var.rand_only, deltaset = false;
    unsigned int n = var.n;
    unsigned int rnum = var.m;
    unsigned int walk_len = var.walk_steps;
    unsigned int n_threads = var.n_threads;
    const NT err = var.err;
    RNGType &rng = var.rng;

    //0. Get the Chebychev ball (largest inscribed ball) with center and radius
    Point c=InnerBall.first;
    NT radius=InnerBall.second;
    if (var.ball_walk){
        if(var.delta<0.0){
            var.delta = 4.0 * radius / NT(n);
            deltaset = true;
        }
    }
    
    //1. Rounding of the polytope if round=true
    NT round_value=1;
    if(round){
        #ifdef VOLESTI_DEBUG
        if(print) std::cout<<"\nRounding.."<<std::endl;
        #endif
        double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::pair<NT,NT> res_round = rounding_min_ellipsoid(P,InnerBall,var);
        round_value=res_round.first;
        double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
        #ifdef VOLESTI_DEBUG
        if(print) std::cout << "Rounding time = " << tstop1 - tstart1 << std::endl;
        #endif
        std::pair<Point,NT> res=P.ComputeInnerBall();
        c=res.first; radius=res.second;
    }

    if (var.ball_walk){
        if(deltaset){
            var.delta = 4.0 * radius / NT(n);
        }
    }

    VT c_e(n);
    for(unsigned int i=0; i<n; i++){
        c_e(i)=c[i];  // write chebychev center in an eigen vector
    }
    P.shift(c_e);
    c=Point(n);

    rnum=rnum/n_threads;
    NT vol=0;
        
        
        
    // Perform the procedure for a number of threads and then take the average
    for(unsigned int t=0; t<n_threads; t++){
        // 2. Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
        #ifdef VOLESTI_DEBUG
        if(print) std::cout<<"\nGenerate the first random point in P"<<std::endl;
        #endif
        
        Point p = get_point_on_Dsphere<RNGType , Point>(n, radius);
        //p=p+c;
        
        std::list<Point> randPoints; //ds for storing rand points
        //use a large walk length e.g. 1000
        
        rand_point_generator(P, p, 1, 50*n, randPoints, var);
        double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
        
        
        // 3. Sample "rnum" points from P
        #ifdef VOLESTI_DEBUG
        if(print) std::cout<<"\nCompute "<<rnum<<" random points in P"<<std::endl;
        #endif
        rand_point_generator(P, p, rnum-1, walk_len, randPoints, var);
        
        double tstop2 = (double)clock()/(double)CLOCKS_PER_SEC;
        #ifdef VOLESTI_DEBUG
        if(print) std::cout << "First random points construction time = " << tstop2 - tstart2 << std::endl;
        #endif

        // 4.  Construct the sequence of balls
        // 4a. compute the radius of the largest ball
        NT current_dist, max_dist=NT(0);
        for(typename  std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
            current_dist=(*pit).squared_length();
            if(current_dist>max_dist){
                max_dist=current_dist;
            }
        }
        max_dist=std::sqrt(max_dist);
        #ifdef VOLESTI_DEBUG
        if(print) std::cout<<"\nFurthest distance from Chebychev point= "<<max_dist<<std::endl;
        #endif

        //
        // 4b. Number of balls
        int nb1 = n * (std::log(radius)/std::log(2.0));
        int nb2 = std::ceil(n * (std::log(max_dist)/std::log(2.0)));

        #ifdef VOLESTI_DEBUG
        if(print) std::cout<<"\nConstructing the sequence of balls"<<std::endl;
        #endif

        std::vector<Ball> balls;
        
        for(int i=nb1; i<=nb2; ++i){

            if(i==nb1){
                balls.push_back(Ball(c,radius*radius));
                vol = (std::pow(M_PI,n/2.0)*(std::pow(balls[0].radius(), n) ) ) / (tgamma(n/2.0+1));
            }else{
                balls.push_back(Ball(c,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
            }

        }
        assert(!balls.empty());

        #ifdef VOLESTI_DEBUG
        if (print) std::cout<<"---------"<<std::endl;
        #endif

        // 5. Estimate Vol(P)

        typename std::vector<Ball>::iterator bit2=balls.end();
        bit2--;

        while(bit2!=balls.begin()){

            //each step starts with some random points in PBLarge stored in list "randPoints"
            //these points have been generated in a previous step

            BallPoly PBLarge(P,*bit2);
            --bit2;
            BallPoly PBSmall(P,*bit2);

            #ifdef VOLESTI_DEBUG
            if(print)
                std::cout<<"("<<balls.end()-bit2<<"/"<<balls.end()-balls.begin()<<") Ball ratio radius="
                        <<PBLarge.second().radius()<<","<<PBSmall.second().radius()<<std::endl;
            #endif

            // choose a point in PBLarge to be used to generate more rand points
            Point p_gen = *randPoints.begin();

            // num of points in PBSmall and PBLarge
            unsigned int nump_PBSmall = 0;
            unsigned int nump_PBLarge = randPoints.size();

            #ifdef VOLESTI_DEBUG
            if(print) std::cout<<"Points in PBLarge="<<randPoints.size()
                              <<std::endl;
            #endif

            //keep the points in randPoints that fall in PBSmall
            typename std::list<Point>::iterator rpit=randPoints.begin();
            while(rpit!=randPoints.end()){
                if (PBSmall.second().is_in(*rpit) == 0){//not in
                    rpit=randPoints.erase(rpit);
                } else {
                    ++nump_PBSmall;
                    ++rpit;
                }
            }

            #ifdef VOLESTI_DEBUG
            if(print) std::cout<<"Points in PBSmall="<<randPoints.size()
                              <<"\nRatio= "<<NT(nump_PBLarge)/NT(nump_PBSmall)
                             <<std::endl;
            #endif

            #ifdef VOLESTI_DEBUG
            if(print) std::cout<<"Generate "<<rnum-nump_PBLarge<<  " more "
                              <<std::endl;
            #endif

            //generate more random points in PBLarge to have "rnum" in total
            rand_point_generator(PBLarge,p_gen,rnum-nump_PBLarge,walk_len,randPoints,PBSmall,nump_PBSmall,var);

            vol *= NT(rnum)/NT(nump_PBSmall);

            #ifdef VOLESTI_DEBUG
            if(print) std::cout<<nump_PBSmall<<"/"<<rnum<<" = "<<NT(rnum)/nump_PBSmall
                              <<"\ncurrent_vol = "<<vol
                            <<"\n--------------------------"<<std::endl;
            #endif

            //don't continue in pairs of balls that are almost inside P, i.e. ratio ~= 2
        }
    }
    vol=round_value*vol;
    #ifdef VOLESTI_DEBUG
    if(print) std::cout<<"rand points = "<<rnum<<std::endl;
    if(print) std::cout<<"walk len = "<<walk_len<<std::endl;
    if(print) std::cout<<"round_value: "<<round_value<<std::endl;
    if(print) std::cout<<"volume computed: "<<vol<<std::endl;
    #endif

    return vol;
}



// Implementation is based on algorithm from paper "A practical volume algorithm",
// Springer-Verlag Berlin Heidelberg and The Mathematical Programming Society 2015
// Ben Cousins, Santosh Vempala
template <class Polytope, class UParameters, class GParameters, class Point, typename NT>
NT volume_gaussian_annealing(Polytope &P,
                             GParameters &var,  // constans for volume
                             UParameters &var2,
                             std::pair<Point,NT> InnerBall) {

    typedef typename Polytope::MT 	MT;
    typedef typename Polytope::VT 	VT;
    typedef typename UParameters::RNGType RNGType;
    const NT maxNT = std::numeric_limits<NT>::max();
    const NT minNT = std::numeric_limits<NT>::lowest();
    NT vol;
    bool round = var.round, done;
    bool print = var.verbose;
    bool rand_only = var.rand_only, deltaset = false;
    unsigned int n = var.n, steps;
    unsigned int walk_len = var.walk_steps, m = P.num_of_hyperplanes();
    unsigned int n_threads = var.n_threads, min_index, max_index, index, min_steps;
    NT error = var.error, curr_eps, min_val, max_val, val, frac = var.frac;
    RNGType &rng = var.rng;
    typedef typename std::vector<NT>::iterator viterator;

    // Consider Chebychev center as an internal point
    Point c=InnerBall.first;
    NT radius=InnerBall.second;
    if (var.ball_walk){
        if(var.delta<0.0){
            var.delta = 4.0 * radius / NT(n);
            var.deltaset = true;
        }
    }

    // rounding of the polytope if round=true
    NT round_value=1;
    if(round){
        #ifdef VOLESTI_DEBUG
        if(print) std::cout<<"\nRounding.."<<std::endl;
        #endif
        double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::pair<NT,NT> res_round = rounding_min_ellipsoid(P,InnerBall,var2);
        double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
        #ifdef VOLESTI_DEBUG
        if(print) std::cout << "Rounding time = " << tstop1 - tstart1 << std::endl;
        #endif
        round_value = res_round.first;
        std::pair<Point,NT> res = P.ComputeInnerBall();
        c = res.first; radius = res.second;
    }

    // Move chebychev center to origin and apply the same shifting to the polytope
    P.shift(Eigen::Map<VT>(&c.get_coeffs()[0], c.dimension()));

    // Scale the polytope s.t. the radius of the chebychev ball equals to 1
    P.set_mat(P.get_mat() * (MT::Identity(n, n) * radius));
    round_value *= std::pow(radius, NT(n));

    // Maybe remove them
    std::pair<Point,NT> res = P.ComputeInnerBall();
    radius = res.second;
    var.che_rad = radius;
    //res.first.print();

    // Initialization for the schedule annealing
    std::vector<NT> a_vals;
    NT ratio = var.ratio, C = var.C;
    unsigned int N = var.N;

    // Computing the sequence of gaussians
    #ifdef VOLESTI_DEBUG
    if(print) std::cout<<"\n\nComputing annealing...\n"<<std::endl;
    #endif
    double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
    get_annealing_schedule(P, radius, ratio, C, frac, N, var, error, a_vals);
    double tstop2 = (double)clock()/(double)CLOCKS_PER_SEC;
    #ifdef VOLESTI_DEBUG
    if(print) std::cout<<"All the variances of schedule_annealing computed in = "<<tstop2-tstart2<<" sec"<<std::endl;
    #endif

    unsigned int mm = a_vals.size()-1, j=0;

    #ifdef VOLESTI_DEBUG
    if(print){
        for (viterator avalIt = a_vals.begin(); avalIt!=a_vals.end(); avalIt++, j++){
            std::cout<<"a_"<<j<<" = "<<*avalIt<<" ";
        }
        std::cout<<"\n"<<std::endl;
    }
    #endif

    // Initialization for the approximation of the ratios
    std::vector<NT> fn(mm,0), its(mm,0), lamdas(m,0);
    unsigned int W = var.W;
    std::vector<NT> last_W2(W,0);
    vol=std::pow(M_PI/a_vals[0], (NT(n))/2.0)*std::abs(round_value);
    Point p(n); // The origin is in the Chebychev center of the Polytope
    Point p_prev=p;
    unsigned int coord_prev, i=0;
    viterator fnIt = fn.begin(), itsIt = its.begin(), avalsIt = a_vals.begin(), minmaxIt;

    #ifdef VOLESTI_DEBUG
    if(print) std::cout<<"volume of the first gaussian = "<<vol<<"\n"<<std::endl;
    if(print) std::cout<<"computing ratios..\n"<<std::endl;
    #endif

    // Compute the first point if CDHR is requested
    if(var.cdhr_walk){
        gaussian_first_coord_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);
    } else if (var.hmc_refl){
        p = get_point_in_Dsphere<RNGType, Point>(P.dimension(), radius);
    }
    for ( ; fnIt != fn.end(); fnIt++, itsIt++, avalsIt++, i++) { //iterate over the number of ratios
        //initialize convergence test
        curr_eps = error/std::sqrt((NT(mm)));
        done=false;
        min_val = minNT;
        max_val = maxNT;
        min_index = W-1;
        max_index = W-1;
        index = 0;
        min_steps=0;
        std::vector<NT> last_W=last_W2;

        // Set the radius for the ball walk if it is requested
        if (var.ball_walk) {
            if (var.deltaset) {
                var.delta = 4.0 * radius / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n));
            }
        }

        while(!done || (*itsIt)<min_steps){

            gaussian_next_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);

            *itsIt = *itsIt + 1.0;
            *fnIt = *fnIt + eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
            val = (*fnIt) / (*itsIt);

            last_W[index] = val;
            if(val<=min_val){
                min_val = val;
                min_index = index;
            }else if(min_index==index){
                minmaxIt = std::min_element(last_W.begin(), last_W.end());
                min_val = *minmaxIt;
                min_index = std::distance(last_W.begin(), minmaxIt);
            }

            if(val>=max_val){
                max_val = val;
                max_index = index;
            }else if(max_index==index){
                minmaxIt = std::max_element(last_W.begin(), last_W.end());
                max_val = *minmaxIt;
                max_index = std::distance(last_W.begin(), minmaxIt);
            }

            if( (max_val-min_val)/max_val<=curr_eps/2.0 ){
                done=true;
            }

            index = index%W+1;

            if(index==W) index=0;
        }
        #ifdef VOLESTI_DEBUG
        if(print) std::cout<<"ratio "<<i<<" = "<<(*fnIt) / (*itsIt)<<" N_"<<i<<" = "<<*itsIt<<std::endl;
        #endif
        vol = vol*((*fnIt) / (*itsIt));
    }
    // Compute and print total number of steps in verbose mode only
    #ifdef VOLESTI_DEBUG
    if (print) {
        NT sum_of_steps = 0.0;
        for(viterator it = its.begin(); it != its.end(); ++it) {
            sum_of_steps += *it;
        }
        steps= int(sum_of_steps);
        std::cout<<"\nTotal number of steps = "<<steps<<"\n"<<std::endl;
    }
    #endif

    return vol;
}


#endif
