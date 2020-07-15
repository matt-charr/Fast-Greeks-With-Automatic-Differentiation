/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#ifndef BLACKSCHOLES_HPP
#define BLACKSCHOLES_HPP

#include "random.hpp"
#include "toolBox.hpp"
#include "features.hpp"

class blackScholes
{
    public:
        blackScholes () = default;
        blackScholes (vec const& s0, vec const& sigma, mat const& corr, double T, double r): 
            s0 (s0), sigma (sigma), mu (r - 0.5 * sigma % sigma), T (T),  Z (corr) {};
        template <typename Tgen>
        features_assets<vec> operator () (Tgen& gen)
        {
            auto z = Z (gen);
            auto s = s0 % exp (mu * T + sqrt (T) * sigma % z);
            auto tgt_r = T * s; 
            auto tgt_sigma = (sqrt (T) * z - T * sigma) % s; 
            auto tgt_s = s / s0; 
            return {s, tgt_s, tgt_sigma, tgt_r};
        }
    private:
        vec s0;
        vec sigma;
        vec mu;
        double T;
        multi_normal Z;
};

class multiBlackScholes
{
    public:
        multiBlackScholes () = default;
        multiBlackScholes (vec const& s0, vec const& sigma, mat const& corr, vec const& dates, double r): 
            n (s0.n_elem), m (dates.n_elem-1), s0 (s0), sigma (sigma), mu (r - 0.5 * sigma % sigma), durations (dates.n_elem-1), Z (corr)
        {
            for (unsigned i (0); i < m; ++i) durations (i) = dates (i+1) - dates (i);
        };
        template <typename Tgen>
        features_assets<mat> operator () (Tgen& gen)
        {
            mat S (n, m), dS_dx (n, m), dS_dsigma (n, m), dS_dr (n, m);
            for (unsigned i (0); i < m; ++i)
            {
                double duration = durations (i);
                auto z = Z (gen);
                auto Ds = exp (mu * duration + sqrt (duration) * sigma % z);
                if (i > 0)
                {
                    dS_dsigma.col (i) = (dS_dsigma.col (i-1) + S.col (i-1) % (sqrt (duration) * z - sigma * duration)) % Ds;
                    dS_dr.col (i) = (dS_dr.col (i-1) + duration * S.col (i-1)) % Ds;
                    S.col (i) = S.col (i-1) % Ds;
                }
                else 
                {
                    auto s0_Ds = s0 % Ds;
                    dS_dsigma.col (0) = (sqrt (duration) * z - sigma * duration) % s0_Ds;
                    dS_dr.col (0) = duration * s0_Ds;
                    S.col (0) = s0_Ds;
                }
            }
            return {S, dS_dx, dS_dsigma, dS_dr};
        }
    private:
        unsigned n, m;
        vec s0;
        vec sigma;
        vec mu;
        vec durations;
        multi_normal Z;
};

#endif
