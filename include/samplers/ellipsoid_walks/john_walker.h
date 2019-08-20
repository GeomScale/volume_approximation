// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2019 Vissarion Fisikopoulos

// Copyright (c) 2018-2019 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef PWALK_JOHN_WALKER_HPP_
#define PWALK_JOHN_WALKER_HPP_

#include <cmath>
#include <Eigen/Dense>
#include "math_functions.h"
#include "walker.h"

//namespace pwalk {

template <typename Dtype>
class JohnWalker: public Walker<Dtype> {
public:
  JohnWalker(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& initialization, const Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>& cons_A, const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& cons_b, const Dtype r) : Walker<Dtype>(initialization, cons_A, cons_b), r_(r), alpha_(1. - 1. / std::log2(2.*Dtype(cons_A.rows())/Dtype(cons_A.cols()))), beta_(Dtype(cons_A.cols())/2./Dtype(cons_A.rows())), curr_weight_(Eigen::Matrix<Dtype, Eigen::Dynamic, 1>::Ones(cons_A.rows())){}

  // getter for radius
  Dtype getRadius() {
    return r_;
  }

  void proposal(Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample){
      Eigen::Matrix<Dtype, Eigen::Dynamic, 1> gaussian_step = Eigen::Matrix<Dtype, Eigen::Dynamic, 1>::Zero(this->nb_dim_);
      sample_gaussian<Dtype>(this->nb_dim_, 0., 1., gaussian_step);

      // get hessian
      Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_sqrt_inv_hess = Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>::Zero(this->nb_dim_, this->nb_dim_);
      sqrtInvHessBarrier(this->curr_sample_, new_sqrt_inv_hess);

      new_sample = this->curr_sample_ + r_ / std::sqrt(Dtype(this->nb_dim_))  * (new_sqrt_inv_hess * gaussian_step);
  }

  bool acceptRejectReverse(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample){
      // get hessian on y
      Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_sqrt_inv_hess_y = Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>::Zero(this->nb_dim_, this->nb_dim_);
      sqrtInvHessBarrier(new_sample, new_sqrt_inv_hess_y);
      // get hessian on x
      Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_sqrt_inv_hess_x = Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>::Zero(this->nb_dim_, this->nb_dim_);
      sqrtInvHessBarrier(this->curr_sample_, new_sqrt_inv_hess_x);

      Dtype scale = r_ / std::sqrt(Dtype(this->nb_dim_));
      Dtype p_y_to_x = gaussian_density<Dtype>(this->curr_sample_, new_sample, new_sqrt_inv_hess_y.inverse()/scale);
      Dtype p_x_to_y = gaussian_density<Dtype>(new_sample, this->curr_sample_, new_sqrt_inv_hess_x.inverse()/scale);

      Dtype ar_ratio = std::min<Dtype>(1., p_y_to_x/p_x_to_y);

      Dtype random_num = rng_uniform<Dtype>(0., 1.);
      // lazy version of the walk
      if (random_num > ar_ratio) {
          return false;
      }

      return true;
  }

  bool doSample(Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample, const Dtype lazy = Dtype(0.5)){
      proposal(new_sample);
      this->nb_curr_samples_ += 1;
      // for lazy markov chain
      Dtype random_num = rng_uniform<Dtype>(0., 1.);
      // check balance and check in polytope
      if (random_num < lazy && this->checkInPolytope(new_sample) && acceptRejectReverse(new_sample)){
          this->curr_sample_ = new_sample;
          return true;
      } else {
          new_sample = this->curr_sample_;
          return false;
      }
  }

  void sqrtInvHessBarrier(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample, Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>& new_sqrt_inv_hess){
      Eigen::Matrix<Dtype, Eigen::Dynamic, 1> inv_slack = (this->cons_b_ - this->cons_A_ * new_sample).cwiseInverse();

      Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> half_hess = inv_slack.asDiagonal()* this->cons_A_;

      Eigen::Matrix<Dtype, Eigen::Dynamic, 1> gradient;
      Eigen::Matrix<Dtype, Eigen::Dynamic, 1> score;
      Eigen::Matrix<Dtype, Eigen::Dynamic, 1> weight_half_alpha;
      Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> weight_half_hess;
      Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_hess;
      Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_hess_inv;
      Eigen::Matrix<Dtype, Eigen::Dynamic, 1> beta_ones = beta_ * Eigen::Matrix<Dtype, Eigen::Dynamic, 1>::Ones(this->nb_cons_);

      Eigen::Matrix<Dtype, Eigen::Dynamic, 1> next_weight = curr_weight_;
      // compute scores using gradient descent
      do {
          curr_weight_ = next_weight;
          weight_half_alpha = curr_weight_.array().pow(alpha_/2.);

          weight_half_hess = (weight_half_alpha.cwiseProduct(inv_slack)).asDiagonal() * this->cons_A_ ;

          new_hess = weight_half_hess.transpose() * weight_half_hess;

          new_hess_inv = new_hess.inverse();

          score = ((weight_half_hess * new_hess_inv).cwiseProduct(weight_half_hess)).rowwise().sum();

          // gradient = Eigen::Matrix<Dtype, Eigen::Dynamic, 1>::Ones(this->nb_cons_) - (score + beta_ones).cwiseQuotient(curr_weight_);

          // curr_weight_ = (curr_weight_ - 0.5 * gradient).cwiseMax(beta_ones);
          next_weight = (0.5*(curr_weight_ + score + beta_ones)).cwiseMax(beta_ones);

      } while ((next_weight - curr_weight_).template lpNorm<Eigen::Infinity>() > Dtype(0.00001));

      // std::cout << "inv_slack" << inv_slack.transpose() << std::endl;
      // std::cout << "score" << score.transpose() << std::endl;
      // std::cout << "curr_weight_ " << curr_weight_.transpose() << std::endl;

      // compute john hessian
      Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> john_new_hess = half_hess.transpose() * curr_weight_.asDiagonal() * half_hess;

      // compute eigenvectors and eigenvalues
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> > es(john_new_hess);

      Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> V = es.eigenvectors();
      Eigen::Matrix<Dtype, Eigen::Dynamic, 1> Dv = es.eigenvalues();
      new_sqrt_inv_hess = V * Dv.cwiseInverse().cwiseSqrt().asDiagonal() * V.transpose();
  }

private:
  const Dtype r_;
  const Dtype alpha_;
  const Dtype beta_;

  Eigen::Matrix<Dtype, Eigen::Dynamic, 1> curr_weight_;
};

//} // namespace pwalk

#endif // PWALK_JOHN_WALKER_HPP_

