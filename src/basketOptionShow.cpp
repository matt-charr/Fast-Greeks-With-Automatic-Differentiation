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

    const unsigned d (16); // Number of assets
    double T = 1. / 3.; // Maturity
    double K = 120; // Strike
    vec w (d); w.fill (1. / d); // Strike

    double r (-0.01); // risk free rate
    vec s0 (d); s0.fill (140); // spots
    vec sigma (d); sigma.fill (0.4); // volatilities
    mat corr (d, d); corr.fill (0.85); for (unsigned i (0); i < d; ++i) corr(i, i) = 1; // Correlations

    auto bs = blackScholes (s0, sigma, corr, T, r);
    auto opt = adj_basketOption (d, w, T, K, r);

    auto mode = 'r'; // d: delta, v: vega, p: price, r: rho
    unsigned MAX = 100000; unsigned batchSize = 100;
    auto X = composed (opt, bs);
    std::ofstream file; file.open ("/Users/MattCharr/Desktop/Algorithmic Differentiation/data/basketOptionShow.txt");
    showMonteCarlo (X, gen, file, mode, MAX, batchSize);

    return 0;
}

