// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef DIAGNOSTICS_ESS_UPDATER_HPP
#define DIAGNOSTICS_ESS_UPDATER_HPP

#include "ess_updater_autocovariance.hpp"


/**
   This is a class that updates the effective sample size (ess) of a sample given a new chain
   using Welford's algorithm to update the average values and the variance estimates where needed.
   The chains has to be of the same length. The ess estimation exploits Geyer's stable estimator
   for the autocovariance and the Geyer's conversion to a monotone sequence, given in,

   Charles J. Geyer, Practical Markov Chain Monte Carlo, Statistical Science 1992.

 * @tparam NT number type
 * @tparam VT vector type
 * @tparam MT matrix type
*/
template <typename NT, typename VT, typename MT>
class ESSestimator {

private:
   unsigned int         num_draws, max_s, s, d, num_chains, jj; 
   VT                   cm_mean, cm_var, cv_mean, draws, var_plus, ess, auto_cov; 
   NT                   oldM, rho_hat_odd, rho_hat_even, mean_var, M2, delta, new_elem;
   MT                   acov_s_mean, rho_hat_s;
   
public:
    ESSestimator() {}

    ESSestimator(unsigned int const& _ndraws, unsigned int const& _dim) 
    {
        num_draws = _ndraws;
        d = _dim;
        num_chains = 0;

        cm_mean.setZero(d);
        cm_var.setZero(d);
        cv_mean.setZero(d);
        var_plus.setZero(d);
        ess.setZero(d);
        draws.setZero(num_draws);
        acov_s_mean.setZero(num_draws-3, d);
        rho_hat_s.setZero(num_draws, d);
   }

    void update_estimator(MT const& samples) 
    {
        num_chains++;
   
        for (int i = 0; i < d; i++)
        {
            draws = samples.row(i).transpose();
            compute_autocovariance<NT>(draws, auto_cov);

            new_elem = draws.mean();
            delta = new_elem - cm_mean.coeff(i);
            cm_mean(i) += delta / NT(num_chains);
            cm_var(i) += delta * (new_elem - cm_mean(i));

            new_elem = auto_cov.coeff(0) * NT(num_draws) / (NT(num_draws) - 1.0);
            delta = new_elem - cv_mean.coeff(i);
            cv_mean(i) += delta / NT(num_chains);

            new_elem = auto_cov.coeff(1);
            delta = new_elem - acov_s_mean.coeff(0, i);
            acov_s_mean(0, i) += delta / NT(num_chains);
            jj = 1;
            while (jj < num_draws-4)
            {
                new_elem = auto_cov.coeff(jj+1);
                delta = new_elem - acov_s_mean.coeff(jj, i);
                acov_s_mean(jj, i) += delta / NT(num_chains);

                new_elem = auto_cov.coeff(jj+2);
                delta = new_elem - acov_s_mean.coeff(jj+1, i);
                acov_s_mean(jj+1, i) += delta / NT(num_chains);

                jj += 2;
            }
        }
    }


    void estimate_effective_sample_size()
    {
        rho_hat_s.setZero(num_draws, d); 

        var_plus = cv_mean * (NT(num_draws) - 1.0) / NT(num_draws);
        if (num_chains > 1) 
        {
            VT cm_var_temp = cm_var * (1.0 / (NT(num_chains)-1.0));
            var_plus += cm_var_temp;
        }

        for (int i = 0; i < d; i++)
        {
            rho_hat_even = 1.0;
            rho_hat_s(0, i) = rho_hat_even;
            rho_hat_odd = 1 - (cv_mean.coeff(i) - acov_s_mean.coeff(0, i)) / var_plus.coeff(i);
            rho_hat_s(1, i) = rho_hat_odd;

            s = 1;
            while (s < (num_draws - 4) && (rho_hat_even + rho_hat_odd) > 0) 
            {
                rho_hat_even = 1.0 - (cv_mean.coeff(i) - acov_s_mean.coeff(s, i)) / var_plus.coeff(i);
                rho_hat_odd = 1.0 - (cv_mean.coeff(i) - acov_s_mean.coeff(s+1, i)) / var_plus.coeff(i);
                if ((rho_hat_even + rho_hat_odd) >= 0) 
                {
                    rho_hat_s(s + 1, i) = rho_hat_even;
                    rho_hat_s(s + 2, i) = rho_hat_odd;
                }
                s += 2;
            }

            max_s = s;
            // this is used in the improved estimate
            if (rho_hat_even > 0) 
            {
                rho_hat_s(max_s + 1, i) = rho_hat_even;
            }

            // Convert Geyer's positive sequence into a monotone sequence
            for (jj = 1; jj <= max_s - 3; jj += 2) 
            {
                if (rho_hat_s(jj + 1, i) + rho_hat_s.coeff(jj + 2, i) > rho_hat_s.coeff(jj - 1, i) + rho_hat_s.coeff(jj, i)) 
                {
                    rho_hat_s(jj + 1, i) = (rho_hat_s.coeff(jj - 1, i) + rho_hat_s.coeff(jj, i)) / 2.0;
                    rho_hat_s(jj + 2, i) = rho_hat_s.coeff(jj + 1, i);
                }
            }
            NT num_total_draws = NT(num_chains) * NT(num_draws);
            NT tau_hat = -1.0 + 2.0 * rho_hat_s.col(i).head(max_s).sum() + rho_hat_s.coeff(max_s + 1, i);
            ess(i) = num_total_draws / tau_hat;
        }
    }


    VT get_effective_sample_size() 
    {
      return ess;
    }

};


#endif

