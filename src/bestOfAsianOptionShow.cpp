/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#include <iostream>

#include "monteCarlo.hpp"
#include "bestOfAsianOption.hpp"
#include "blackScholes.hpp"

#include "toolBox.hpp"

int main ()
{
    init_seed ();
    auto gen = genrand_real3;

    const unsigned n (3);
    const unsigned m (4);

    double T = 1;
    double K = 0.1;
    auto dates = arma::linspace (0, T, m);

    double r (-0.01);
    vec s0 (n); s0.fill (1);
    vec sigma (n); sigma.fill (0.3);
    mat corr (n, n); corr.fill (0.5); for (unsigned i (0); i < n; ++i) corr(i, i) = 1;

    auto bs = multiBlackScholes (s0, sigma, corr, dates, r);
    auto opt = tgt_bestOfAsianOption (n, m, T, K, r);

    auto mode = 'r';
    unsigned MAX = 100000;
    unsigned batchSize = 100;
    auto X = composed (opt, bs);
    std::ofstream file; file.open ("/Users/MattCharr/Desktop/Algorithmic Differentiation/data/bestOfAsianOptionShow.txt");
    showMonteCarlo (X, gen, file, mode, MAX, batchSize);

    return 0;
}

