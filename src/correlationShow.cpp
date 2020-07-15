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

    const unsigned n (2);

    double T = 1;
    double K = 100;
    vec w (n); w.fill (1. / n);

    vec s0 (n); s0.fill (100);
    vec sigma (n); sigma.fill (0.3);
    mat corr (n, n); corr.fill (0.5); for (unsigned i (0); i < n; ++i) corr(i, i) = 1;
    
    std::ofstream file; file.open ("/Users/MattCharr/Desktop/Algorithmic Differentiation/data/basketOptionCorrelationShow.txt");
    auto opt = basketOptionCorrelation (n, s0, sigma, corr, T, K, w);
    auto mv = opt.monteCarloShow<basketOptionCorrelation> (gen, file, 1000, 1000000, 0, 1);

    std::cout << mv << std::endl;

    return 0;
}

