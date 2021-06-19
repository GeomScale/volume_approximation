// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ODE_SOLVERS_LEAPFROG_HPP
#define ODE_SOLVERS_LEAPFROG_HPP

template <
typename Point,
typename NT,
typename Polytope,
typename func
>
struct LeapfrogODESolver {

  struct update_parameters
{
        update_parameters()
                :   facet_prev(0), hit_ball(false), inner_vi_ak(0.0), ball_inner_norm(0.0)
        {}
        int facet_prev;
        bool hit_ball;
        double inner_vi_ak;
        double ball_inner_norm;
};

  update_parameters _update_parameters;

  typedef std::vector<Point> pts;

  typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

  typedef std::vector<Polytope*> bounds;
  typedef typename Polytope::VT VT;

  unsigned int dim;

  std::vector<VT> Ar, Av;
  std::vector<NT> lambda_prev;

  NT eta;
  NT eta0;
  NT t;
  NT dl = 0.995;

  func F;
  bounds Ks;

  // Contains the sub-states
  pts xs;
  pts xs_prev;

  MT _AA;

  std::pair<NT, int> pbpair;

  unsigned long long num_reflections = 0;
  unsigned long long num_steps = 0;

  bool adaptive = true;

  LeapfrogODESolver(NT initial_time, NT step, pts initial_state, func oracle, bounds boundaries, bool adaptive_=true) :
  eta(step), eta0(step), t(initial_time), F(oracle), Ks(boundaries), xs(initial_state), adaptive(adaptive_) {
    dim = xs[0].dimension();
    _update_parameters = update_parameters();
    initialize();
  };



  void initialize() {
      for (unsigned int i = 0; i < xs.size(); i++) {
          VT ar, av;
          if (Ks[i] != NULL) {
              ar.setZero(Ks[i]->num_of_hyperplanes());
              ar = Ks[i]->get_mat() * xs[i].getCoefficients();
              av.setZero(Ks[i]->num_of_hyperplanes());
              _AA.noalias() = Ks[i]->get_mat() * Ks[i]->get_mat().transpose();
              Ks[i]->resetFlags();
          }
          Ar.push_back(ar);
          Av.push_back(av);
          lambda_prev.push_back(NT(0));
      }
      //step();
  }

  void disable_adaptive() {
      adaptive = false;
  }

  void enable_adaptive() {
      adaptive = true;
  }

  void step(int k, bool accepted) {
    num_steps++;

    if (adaptive) eta = (eta0 * num_steps) / (num_steps + num_reflections);

    xs_prev = xs;
    unsigned int x_index, v_index, it;
    t += eta;
    for (unsigned int i = 1; i < xs.size(); i += 2) {
        //pbpair.second  = -1;
      x_index = i - 1;
      v_index = i;

      // v' <- v + eta / 2 F(x)
      Point z = F(v_index, xs_prev, t);
      z = (eta / 2) * z;
      xs[v_index] = xs[v_index] + z;

      // x <- x + eta v'
      Point y = xs[v_index];

      if (Ks[x_index] == NULL) {
        xs[x_index] = xs_prev[x_index] + eta*y;
      }
      else {
        // Find intersection (assuming a line trajectory) between x and y
        bool step_completed = false;
        NT T = eta;
        if (k == 0 && !accepted) {
          Ar[x_index] = Ks[x_index]->get_mat() * xs_prev[x_index].getCoefficients();
          lambda_prev[x_index] = 0.0;
          Ks[x_index]->resetFlags();
        }

        pbpair =  Ks[x_index]->line_positive_intersect(xs_prev[x_index], y, Ar[x_index], Av[x_index],
                                                       lambda_prev[x_index], _update_parameters);
        if (T <= pbpair.first) {
          xs[x_index] = xs_prev[x_index] + T * y;
          xs[v_index] = y;
          lambda_prev[x_index] = T;
          Ks[x_index]->update_position_internal(T);
          step_completed = true;
        }

        if (!step_completed) {
          lambda_prev[x_index] = dl * pbpair.first;
          xs_prev[x_index] = xs_prev[x_index] + lambda_prev[x_index] * y;
          T -= lambda_prev[x_index];
          Ks[x_index]->update_position_internal(lambda_prev[x_index]);
          Ks[x_index]->compute_reflection(y, xs_prev[x_index], _update_parameters);
          num_reflections++;

          while (true)
          {
            pbpair =  Ks[x_index]->line_positive_intersect(xs_prev[x_index], y, Ar[x_index], Av[x_index],
                                                           lambda_prev[x_index], _AA, _update_parameters);
            if (T <= pbpair.first) {
              xs[x_index] = xs_prev[x_index] + T * y;
              xs[v_index] = y;
              lambda_prev[x_index] = T;
              Ks[x_index]->update_position_internal(T);
              break;
            }
            lambda_prev[x_index] = dl * pbpair.first;
            xs_prev[x_index] = xs_prev[x_index] + lambda_prev[x_index] * y;
            T -= lambda_prev[x_index];
            Ks[x_index]->update_position_internal(lambda_prev[x_index]);
            Ks[x_index]->compute_reflection(y, xs_prev[x_index], _update_parameters);
            num_reflections++;
          }
        }
      }

      // tilde v <- v + eta / 2 F(tilde x)
      z = F(v_index, xs, t);
      z = (eta / 2) * z;
      xs[v_index] = xs[v_index] + z;

    }

  }

  void print_state() {
    for (int j = 0; j < xs.size(); j ++) {
      for (unsigned int i = 0; i < xs[j].dimension(); i++) {
        std::cout << xs[j][i] << " ";
      }
    }
    std::cout << std::endl;
  }

  void steps(int num_steps, bool accepted) {
    for (int i = 0; i < num_steps; i++) step(i, accepted);
  }

  Point get_state(int index) {
    return xs[index];
  }

  void set_state(int index, Point p) {
    xs[index] = p;
  }

};


#endif
