/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#include <iostream>

#include "tgt_double.hpp"

tgt_double f (tgt_double const& x, tgt_double const& y, tgt_double const& z) 
{
    return (log (x * y) + cos (x) * sin (y) * sin (z)) / exp (x);
}

int main ()
{
    double x, y, z;
    double x_dot, y_dot, z_dot;

    std::cout << "x: "; std::cin >> x;
    std::cout << "x coefficient : "; std::cin >> x_dot;
    std::cout << "y: "; std::cin >> y;
    std::cout << "y coefficient : "; std::cin >> y_dot;
    std::cout << "z: "; std::cin >> z;
    std::cout << "z coefficient : "; std::cin >> z_dot; 
    std::cout << std::endl;

    auto tgt_x = make0<tgt_double> (x, {x_dot});
    auto tgt_y = make0<tgt_double> (y, {y_dot});
    auto tgt_z = make0<tgt_double> (z, {z_dot});
    auto tgt = f (tgt_x, tgt_y, tgt_z);

    std::cout << "f (" << x << "," << y << "," << z << ") = " << tgt.get_value () << std::endl;
    std::cout << x_dot << "df_dx + " << y_dot << "df_dy + " << z_dot << "df_dz =" << tgt.get_derivative () << std::endl;

    return 0;
}

