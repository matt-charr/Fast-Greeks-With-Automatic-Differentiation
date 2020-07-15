/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#ifndef FEATURES_HPP
#define FEATURES_HPP

#include "monteCarlo.hpp"

#include "toolBox.hpp"

template<typename Type>
struct features_assets
{
    Type s;
    Type tgt_s;
    Type tgt_sigma;
    Type tgt_r;
};

struct features_option
{
    double p;
    double dp_dr;
    vec dp_ds;
    vec dp_dsigma;
    friend features_option operator% (features_option const& u, features_option const& v)
    {
        return {u.p * v.p, u.dp_dr * v.dp_dr, u.dp_ds % v.dp_ds, u.dp_dsigma % v.dp_dsigma};
    }
    friend features_option operator- (features_option const& u, features_option const& v)
    {
        return {u.p - v.p, u.dp_dr - v.dp_dr, u.dp_ds - v.dp_ds, u.dp_dsigma - v.dp_dsigma};
    }
    friend features_option sqrt (features_option const& u)
    {
        return {sqrt (u.p), sqrt (u.dp_dr), sqrt (u.dp_ds), sqrt (u.dp_dsigma)};
    }
    friend features_option operator* (double k, features_option const& u)
    {
        return {k * u.p, k * u.dp_dr, k * u.dp_ds, k * u.dp_dsigma};
    }
    friend std::ostream& operator<< (std::ostream& os, features_option const& x)
    {
        os << "Price: " << x.p << std::endl;
        os << "Deltas: " << std::endl << x.dp_ds;
        os << "Vegas: " << std::endl << x.dp_dsigma;
        os << "Rho: " << x.dp_dr << std::endl;
        return os;
    }
    friend features_option operator/ (features_option const& u, double n)
    {
        return {u.p / n, u.dp_dr / n, u.dp_ds / n, u.dp_dsigma / n};
    }
    features_option operator+= (features_option const& other)
    {
        p += other.p;
        dp_dr += other.dp_dr;
        dp_ds += other.dp_ds;
        dp_dsigma += other.dp_dsigma;
        return *this;
    }
    features_option operator%= (features_option const& other)
    {
        p *= other.p;
        dp_dr *= other.dp_dr;
        dp_ds %= other.dp_ds;
        dp_dsigma %= other.dp_dsigma;
        return *this;
    }
};

template <typename Tdistribution, typename Tgen>
void showMonteCarlo (Tdistribution& X, Tgen& gen, std::ofstream& os, char mode, unsigned MAX, unsigned batchSize)
{
    mean_var<features_option> mv = monteCarlo<features_option, Tdistribution, Tgen>  (X, batchSize, gen);
    unsigned cpt (batchSize);
    while (cpt < MAX)
    {
        mv += monteCarlo<features_option, Tdistribution, Tgen>  (X, batchSize, gen);
        cpt += batchSize;
        auto mean = mv.mean ();
        auto icSize = mv.icSize ();
        double mean_ (0), icSize_ (0);
        switch (mode)
        {
            case 'p':
                mean_ = mean.p;
                icSize_ = icSize.p;
            break;
            case 'r':
                mean_ = mean.dp_dr;
                icSize_ = icSize.dp_dr;
            break;
            case 'd':
                mean_ = mean.dp_ds (0);
                icSize_ = icSize.dp_ds (0);
            break;
            case 'v':
                mean_ = mean.dp_dsigma (0);
                icSize_ = icSize.dp_dsigma (0);
            break;
        } 
        os << cpt << ";" << mean_ << ";" << icSize_ << std::endl;
    }
    std::cout << mv << std::endl;
}

#endif
