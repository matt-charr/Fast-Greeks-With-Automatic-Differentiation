/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#include <iostream>

#include "monteCarlo.hpp"
#include "basketOption.hpp"
#include "blackScholes.hpp"

#include "toolBox.hpp"

int main ()
{
    init_seed ();
    auto gen = genrand_real3;

    unsigned minimum_number_assets (100);
    unsigned maximum_number_assets (500);
    unsigned batchSize (100);

    double T = 1;
    double K = 100;

    for (unsigned d (minimum_number_assets); d <= maximum_number_assets; d += batchSize)
    {
        vec w (d); w.fill (1. / d);

        double r (-0.01);
        vec s0 (d); s0.fill (110);
        vec sigma (d); sigma.fill (0.3);
        mat corr (d, d); corr.fill (0.5); for (unsigned i (0); i < d; ++i) corr(i, i) = 1;

        auto bs = blackScholes (s0, sigma, corr, T, r);
        
        auto fdm = std_basketOption (d, w, T, K, r);
        auto tm = tgt_basketOption (d, w, T, K, r);
        auto am = adj_basketOption (d, w, T, K, r);
        
        auto std = composed (fdm, bs);
        auto tgt = composed (tm, bs);
        auto adj = composed (am, bs);

        unsigned sampleSize (1000);
        
        auto std_time = timeMonteCarlo<features_option> (std, sampleSize, gen);
        auto tgt_time = timeMonteCarlo<features_option> (tgt, sampleSize, gen);
        auto adj_time = timeMonteCarlo<features_option> (adj, sampleSize, gen);

        std::cout << "Number of assets: " << d << std::endl;
        std::cout << "Time Finite Difference Method: " << std_time << std::endl;
        std::cout << "Time Tangent Mode: " << tgt_time << std::endl;
        std::cout << "Time Adjoint Mode: " << adj_time << std::endl;
    }
    
    return 0;
}

