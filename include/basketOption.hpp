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
class basketOption
{
    public:
        basketOption (unsigned d, vec const& w, double T, double K, double r):
            d (d), payoff ({T, K, w}), r (r), phi (0), dphi_dr (0), dphi_ds (d) {}; 
        virtual features_option operator () (features_assets<vec> const& assets) = 0;
    protected:
        double phi; double dphi_dr; vec dphi_ds;
        unsigned d;
        double r;
        struct
        {
            double T;
            double K;
            vec w;
            inline type_double operator () (type_double const& r, vector<type_double> const& S) 
            { 
                return exp (-r * T) * max (dot (w, S) - K, 0); 
            } 
        } payoff;
};

class std_basketOption : public basketOption<std_double>
{
    public: 
        std_basketOption (unsigned d, vec const& w, double T, double K, double r): 
            basketOption<std_double> (d, w, T, K, r), epsilon (1e-2), vec_epsilon (d, arma::fill::zeros) {};
        features_option operator () (features_assets<vec> const& assets)
        {
            auto R = make0<std_double> (r);
            auto S = make1<std_double> (assets.s);
            phi = payoff (R, S);
            for (unsigned i (0); i < d; ++i) 
            {
                vec_epsilon (i) = epsilon;
                auto S_ = make1<std_double> (assets.s + vec_epsilon);
                dphi_ds (i) = (payoff (R, S_) - phi) / epsilon;
                vec_epsilon (i) = 0;
            }
            dphi_dr = (payoff (R + epsilon, S) - phi) / epsilon;

            auto p = phi;
            auto dp_ds = assets.tgt_s % dphi_ds;
            auto dp_dsigma = assets.tgt_sigma % dphi_ds;
            auto dp_dr = dot (assets.tgt_r, dphi_ds) + dphi_dr;

            return {p, dp_dr, dp_ds, dp_dsigma};
        }
    private:
        double epsilon; vec vec_epsilon;
};

class tgt_basketOption: public basketOption<tgt_double>
{
    public: 
        tgt_basketOption (unsigned d, vec const& w, double T, double K, double r): 
            basketOption<tgt_double> (d, w, T, K, r), r_dot (2*d+1), s_dot (2*d+1, d, arma::fill::zeros), twice_d (2*d) 
            { 
                for (unsigned i (0); i < twice_d; ++i) r_dot (i) = 0;
                r_dot (twice_d) = 1; 
            };        
        features_option operator () (features_assets<vec> const& assets)
        {
            
            for (unsigned i (0); i < d; ++i) s_dot (i, i) = assets.tgt_s (i);
            for (unsigned i (0); i < d; ++i) s_dot (i+d, i) = assets.tgt_sigma (i);
            for (unsigned j (0); j < d; ++j) s_dot (twice_d, j) = assets.tgt_r (j);

            auto R = make0<tgt_double> (r, r_dot);
            auto S = make1<tgt_double> (assets.s, s_dot);
            auto P = payoff (R, S);
            auto D = P.get_derivative ();

            auto p = P.get_value ();
            auto dp_ds = D.head (d);
            auto dp_dsigma = D.subvec (d, twice_d-1);
            auto dp_dr = D (twice_d);

            return {p, dp_dr, dp_ds, dp_dsigma};
        }
    private:
        vec r_dot; mat s_dot;
        unsigned twice_d;
};

class adj_basketOption : public basketOption<adj_double>
{
    public: 
        adj_basketOption (unsigned d, vec const& w, double T, double K, double r): 
            basketOption<adj_double> (d, w, T, K, r) {}
        features_option operator () (features_assets<vec> const& assets)
        { 
            auto t = new Tape ();
            auto R = t->make0<adj_double> (r);
            auto S = t->make1<adj_double> (assets.s);
            auto P = payoff (R, S);
            auto D = P.grad ();
            phi = P.get_value ();
            for (unsigned i (0); i < d; ++i) dphi_ds (i) = D.wrt (S (i));
            dphi_dr = D.wrt (R);

            auto p = phi;
            auto dp_ds = assets.tgt_s % dphi_ds;
            auto dp_dsigma = assets.tgt_sigma % dphi_ds;
            auto dp_dr = dot (assets.tgt_r, dphi_ds) + dphi_dr;

            return {p, dp_dr, dp_ds, dp_dsigma};
        }
};

#endif
