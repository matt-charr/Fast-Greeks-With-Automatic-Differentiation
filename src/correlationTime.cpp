/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#include <iostream>

#include "correlation.hpp"

int main ()
{
    init_seed ();
    auto gen = genrand_real3;
    
    double T = 1;
    double K = 100;

    for (unsigned d (1); d <= 50; d += 10)
    {    
        vec w (d); w.fill (1. / d);
        vec s0 (d); s0.fill (110);
        vec sigma (d); sigma.fill (0.3);
        mat corr (d, d); corr.fill (0.5); for (unsigned i (0); i < d; ++i) corr(i, i) = 1;

        auto opt = basketOptionCorrelation (d, s0, sigma, corr, T, K, w);

        auto time = opt.monteCarloTime (1000, gen);

        std::cout << "Number of assets: " << d << std::endl;
        std::cout << "Time Adjoint Mode: " << time << std::endl;
    }
    
    return 0;
}

