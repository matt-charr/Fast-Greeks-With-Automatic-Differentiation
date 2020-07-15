/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#include <iostream>

#include "adj_double.hpp"

adj_double f (adj_double const& x, adj_double const& y, adj_double const& z) 
{
    return x * y;
}

int main ()
{
    double x, y, z, seed;

    std::cout << "x: "; std::cin >> x;
    std::cout << "y: "; std::cin >> y;
    std::cout << "z: "; std::cin >> z;
    std::cout << "seed: "; std::cin >> seed; std::cout << std::endl;

    auto t = new Tape ();
    auto adj_x = t->make0<adj_double> (x);
    auto adj_y = t->make0<adj_double> (y);
    auto adj_z = t->make0<adj_double> (z);
    auto adj = f (adj_x, adj_y, adj_z);
    auto grad = adj.grad (seed);

    std::cout << "f (" << x << "," << y << "," << z << ") = " << adj.get_value () << std::endl;
    std::cout << "df_dx: (scaled by seed) " << grad.wrt (adj_x) << std::endl;
    std::cout << "df_dy: (scaled by seed) " << grad.wrt (adj_y) << std::endl;
    std::cout << "df_dz: (scaled by seed) " << grad.wrt (adj_z) << std::endl;

    return 0;
}

