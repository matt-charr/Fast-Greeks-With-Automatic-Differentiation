/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#ifndef STDDOUBLE_HPP
#define STDDOUBLE_HPP

#include "toolBox.hpp"

using std_double = double;

// Compute the maximum between x and c
std_double max (std_double x, double c) 
{ 
    if (x > c) return x; else return c; 
}

// Compute the maximum of the elements of the vector x
unsigned max (vec const& x)
{
    unsigned M (0);
    unsigned cpt (0);
    for (auto& x0: x)
    {
        if (x0 > x (M)) M = cpt;
        cpt++;
    } 
    return cpt;
}

// Compute the inner product between w and u
std_double dot (vec const& w, vector<std_double> const& u)
{
    unsigned d (w.n_elem);
    auto sum = w (0) * u (0); 
    for (unsigned i (1); i < d; ++i) sum += w (i) * u (i);
    return sum;
}

// Compute the maximum of the elemtns of the vector u
std_double max (vector<std_double> const& u)
{
    unsigned d (u.n_elem);
    double max (-1e6), val (0);
    for (unsigned i (0); i < d; ++i)
    {
        val = u (i);
        if (val > max) max = val;
    }
    return max;
}

// Compute the vector of the maximum of the column of the matrix U
vector<std_double> max (matrix<std_double> const& U)
{
    unsigned d (U.n_elem);
    vector<std_double> u (d);
    for (unsigned i (0); i < d; ++i) u (i) = max (U (i));
    return u;
}

// Compute the mean of the vector u
std_double mean (vector<std_double> const& u)
{
    unsigned d (u.n_elem);
    double sum (0);
    for (unsigned i (0); i < d; ++i) sum += u (i);
    return sum / d;
}

// Build a std_double
template<typename std_double>
std_double make0 (double s)
{
    return s;
}

// Build a vector of std_double
template<typename std_double>
vector<std_double> make1 (vec const& s)
{
    unsigned d (s.n_elem);
    vector<std_double> res (d);
    for (unsigned i (0); i < d; ++i) res (i) = s (i);
    return res;
}

// Build a matrix of std_double
template<typename std_double>
matrix<std_double> make2 (mat const& s)
{
    unsigned m (s.n_cols);
    matrix<std_double> res (m);
    for (unsigned j (0); j < m; ++j) res (j) = make1<std_double> (s.col (j));
    return res;
}

#endif
