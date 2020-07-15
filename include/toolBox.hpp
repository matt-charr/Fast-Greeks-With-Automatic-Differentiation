/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#ifndef TOOLBOX_HPP
#define TOOLBOX_HPP

#include "armadillo.hpp"

using uvec = arma::uvec;
using vec = arma::vec;
using sp_mat = arma::sp_mat;
using mat = arma::mat;
using cube = arma::cube;

template<typename type_double>
using vector = arma::field<type_double>;

template<typename type_double>
using matrix = arma::field<vector<type_double>>;

#endif
