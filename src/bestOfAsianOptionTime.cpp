/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#include <iostream>

#include "blackScholes.hpp"
#include "monteCarlo.hpp"
#include "bestOfAsianOption.hpp"

#include "toolBox.hpp"

int main ()
{
    init_seed ();
    auto gen = genrand_real3;

    const unsigned m (12);
    double T = 1;
    double K = 0.5;
    auto dates = arma::linspace (0, T, m);
    double r (0.1);

    unsigned minimum_number_assets (100);
    unsigned maximum_number_assets (500);
    unsigned batchSize (100);

    for (unsigned n (minimum_number_assets); n <= maximum_number_assets; n += batchSize)
    {
        vec s0 (n, arma::fill::ones);
        vec sigma (n); sigma.fill (0.3);
        mat corr (n, n); corr.fill (0.5); for (unsigned i (0); i < n; ++i) corr(i, i) = 1;

        auto bs = multiBlackScholes (s0, sigma, corr, dates, r);
    
        auto fdm = std_bestOfAsianOption (n, m, T, K, r);
        auto tm = tgt_bestOfAsianOption (n, m, T, K, r);
        auto am = adj_bestOfAsianOption (n, m, T, K, r);

        auto std = composed (fdm, bs);
        auto tgt = composed (tm, bs);
        auto adj = composed (am, bs);

        unsigned sampleSize (100);

        auto std_time = timeMonteCarlo<features_option> (std, sampleSize, gen);
        auto tgt_time = timeMonteCarlo<features_option> (tgt, sampleSize, gen);
        auto adj_time = timeMonteCarlo<features_option> (adj, sampleSize, gen);

        std::cout << "Number of assets: " << n << std::endl;
        std::cout << "Time Finite Difference Method: " << std_time << std::endl;
        std::cout << "Time Tangent Mode: " << tgt_time << std::endl;
        std::cout << "Time Adjoint Mode: " << adj_time << std::endl;
    }

    return 0;
}

