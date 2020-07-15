/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#ifndef BASKETOPTION_HPP
#define BASKETOPTION_HPP

#include "std_double.hpp"
#include "tgt_double.hpp"
#include "adj_double.hpp"

#include "monteCarlo.hpp"

#include "blackScholes.hpp"

#include "features.hpp"

template<typename type_double>
class bestOfAsianOption
{
    public:
        bestOfAsianOption (unsigned n, unsigned m, double T, double K, double r): phi (0), dphi_dr (0), dphi_ds (n, m-1),
            n (n), m (m-1), r (r), payoff ({T, K}) {}; 
        virtual features_option operator () (features_assets<mat> const& assets) = 0;
    protected:
        double phi; double dphi_dr; mat dphi_ds;
        unsigned n; unsigned m;
        double r;
        struct
        {
            double T;
            double K;
            inline type_double operator () (type_double const& r, matrix<type_double> const& S) 
            { 
                return exp (-r * T) * max (mean (max (S)) - K, 0); 
            } 
        } payoff;
};

class std_bestOfAsianOption: public bestOfAsianOption<std_double>
{
    public: 
        std_bestOfAsianOption (unsigned n, unsigned m, double T, double K, double r):
            bestOfAsianOption<std_double> (n, m, T, K, r), epsilon (1e-2), mat_epsilon (n, m-1, arma::fill::zeros) {};
        features_option operator () (features_assets<mat> const& assets)
        {
            // Compute value and derivatives
            auto R = make0<std_double> (r);
            auto S = make2<std_double> (assets.s);
            phi = payoff (R, S);
            for (unsigned j (0); j < m; ++j)
                for (unsigned i (0); i < n; ++i) 
                {
                    mat_epsilon (i, j) = epsilon;
                    auto S_ = make2<std_double> (assets.s + mat_epsilon);
                    dphi_ds (i, j) = (payoff (R, S_) - phi) / epsilon;
                    mat_epsilon (i, j) = 0;
            }
            dphi_dr = (payoff (R + epsilon, S) - phi) / epsilon;

            // Compute price and sensitivities
            auto p = phi;
            auto dp_ds = arma::zeros (n);
            auto dp_dr = dot (assets.tgt_r, dphi_ds) + dphi_dr;
            vec dp_dsigma (n);
            for (unsigned i (0); i < n; ++i) 
                dp_dsigma (i) = dot (assets.tgt_sigma.row (i), dphi_ds.row (i));
    
            return {p, dp_dr, dp_ds, dp_dsigma};
        }
    private:
        double epsilon; mat mat_epsilon;
};

class tgt_bestOfAsianOption: public bestOfAsianOption<tgt_double>
{
    public: 
        tgt_bestOfAsianOption (unsigned n, unsigned m, double T, double K, double r):
            bestOfAsianOption<tgt_double> (n, m, T, K, r), r_dot (n+1), s_dot (n+1, n, m)
        {
            for (unsigned i (0); i < n; ++i) r_dot (i) = 0;
            r_dot (n) = 1; 
        }
        features_option operator () (features_assets<mat> const& assets)
        {
            for (unsigned k (0); k < m; ++k)
            {
                for (unsigned j (0); j < n; ++j)
                {
                    for (unsigned i (0); i < n; ++i) 
                        if (i == j) s_dot (i, j, k) = assets.tgt_sigma (j, k); else s_dot (i, j, k) = 0;
                    s_dot (n, j, k) = assets.tgt_r (j, k);
                }
            }
            
            // Compute value and derivatives
            auto R = make0<tgt_double> (r, r_dot);
            auto S = make2<tgt_double> (assets.s, s_dot);
            auto P = payoff (R, S);
            auto D = P.get_derivative ();

            // Compute price and sensitivites
            auto p = P.get_value ();
            auto dp_ds = arma::zeros (n);
            auto dp_dsigma = D.head (n);
            auto dp_dr = D (n);

            return {p, dp_dr, dp_ds, dp_dsigma};
        }
    private:
        vec r_dot; cube s_dot;
};

class adj_bestOfAsianOption: public bestOfAsianOption<adj_double>
{
    public: 
        adj_bestOfAsianOption (unsigned n, unsigned m, double T, double K, double r):
            bestOfAsianOption<adj_double> (n, m, T, K, r) {}
        features_option operator () (features_assets<mat> const& assets)
        {
            // Compute value and derivatives
            auto t = new Tape ();
            auto R = t->make0<adj_double> (r);
            auto S = t->make2<adj_double> (assets.s);
            auto P = payoff (R, S);
            auto D = P.grad ();
            phi = P.get_value ();
            for (unsigned j (0); j < m; ++j)
                for (unsigned i (0); i < n; ++i)
                    dphi_ds (i, j) = D.wrt (S (j) (i));
            dphi_dr = D.wrt (R);
            
            // Compute price and sensitivities
            auto p = phi;
            auto dp_ds = arma::zeros (n);
            auto dp_dr = dot (assets.tgt_r, dphi_ds) + dphi_dr;
            vec dp_dsigma (n);
            for (unsigned i (0); i < n; ++i) 
                dp_dsigma (i) = dot (assets.tgt_sigma.row (i), dphi_ds.row (i));
    
            return {p, dp_dr, dp_ds, dp_dsigma};
        }
};

#endif
