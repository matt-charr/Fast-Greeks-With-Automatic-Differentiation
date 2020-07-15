/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#ifndef CORRELATION_HPP
#define CORRELATION_HPP

#include "toolBox.hpp"
#include "random.hpp"
#include "adj_double.hpp"
#include "monteCarlo.hpp"

class basketOptionCorrelation
{
    public:
        basketOptionCorrelation () = default;
        basketOptionCorrelation (unsigned d, vec const& s0, vec const& sigma, mat const& rho, double T, double K, vec w): 
            T (T), d (d), s0 (s0), sigma (sigma), mu (-0.5 * sigma % sigma), Z (d), t (new Tape ()), payoff ({K, w})
        {
            adj_rho_L = t->make_correlation_matrix<adj_double> (rho);
            adj_cholesky (adj_rho_L, t);
        };
        template <typename Tgen>
        sp_mat operator () (Tgen& gen)
        {
            auto cpy = new Tape (t);
            for (auto& x_j: adj_rho_L) for (auto& x_i: x_j) x_i.set_tape (cpy);
            auto z = Z (gen);
            auto z_corr = adj_rho_L * z;
            auto w = sqrt (T) * z_corr;
            auto s = s0 % exp (sigma % w + mu * T);
            auto p = payoff (s);
            auto grad = p.grad ();
            sp_mat dp_drho (d, d);
            for (unsigned j (1); j < d; ++j)
                for (unsigned i (0); i < j; ++i)
                    dp_drho (i, j) = grad.wrt (adj_rho_L (j) (i));
            dp_drho += dp_drho;
            return dp_drho;
        }
        template <typename Tdistribution, typename Tgen>
        mean_var<sp_mat> monteCarloShow (Tgen& gen, std::ofstream& os, unsigned batch_size, unsigned max_size, unsigned i, unsigned j)
        {
            mean_var<sp_mat> mv = monteCarlo<sp_mat, Tdistribution, Tgen>  (*this, batch_size, gen);
            unsigned cpt (batch_size);
            while (cpt < max_size)
            {
                mv += monteCarlo<sp_mat, Tdistribution, Tgen>  (*this, batch_size, gen);
                cpt += batch_size;
                auto mean = mv.mean ();
                auto icSize = mv.icSize ();
                os << cpt << ";" << mean (i, j) << ";" << icSize (i, j) << std::endl;
            }
            return mv;
        }   
        template <typename Tgen>
        float monteCarloTime (unsigned n, Tgen& gen)
        {
            return timeMonteCarlo<sp_mat> (*this, n, gen);
        }
    private:
        unsigned d;
        double T;
        vec s0, sigma, mu;
        multi_normal_iid Z;
        Tape* t;
        matrix<adj_double> adj_rho_L;
        struct 
        {
            double K;
            vec w;
            adj_double operator () (vector<adj_double> const& s)
            {
                return max (dot (w, s) - K, 0);
            }
        } payoff;
};

#endif
