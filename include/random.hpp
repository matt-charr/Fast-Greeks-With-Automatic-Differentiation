/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <cmath>

#include "mt19937ar.hpp"

#include "toolBox.hpp"

void init_seed (unsigned int seed = (unsigned int) std::time (NULL)) 
{
    init_genrand (seed); 
};

template<typename T>
class random_variable
{
    protected:
        T value;
    public:
        using result_type = T;
        random_variable () = default;
        random_variable (T const& value) : value (value) {}
        virtual ~random_variable () {};
        result_type current () const  { return (value); };
};

class uniform: public random_variable<double>
{
    private:
        double left, size;
    public:    
        uniform (double left = 0, double right = 1) : left(left), size(right - left) {};
        template<typename Tgen>
        double operator() (Tgen& gen) { return value = left + size * gen (); };
};

class normal: public random_variable<double>
{
    private:
        double mean, std, U, V, sq_R, scale;
        uniform unif;
        bool flag;
    public:
        normal (double mean = 0, double std = 1): mean (mean), std (std), unif (-1, 1), flag (true) {};
        template<typename Tgen>
        double operator () (Tgen& gen)
        {
            flag = !flag;
            if (!flag)
            {
                do
                {
                    U = unif (gen);
                    V = unif (gen);
                    sq_R = U * U + V * V;
                } while (sq_R > 1.);
                scale = sqrt (-2. * log (sq_R) / sq_R);
                return value = mean + std * U * scale;
            }
            else return value = mean + std * V * scale;
        };
};

class multi_normal_iid: public random_variable<vec>
{
    private:
        unsigned d;
        normal Z;
    public:
        multi_normal_iid (unsigned d): Z (0, 1), d (d) {};
        template<typename Tgen>
        vec operator () (Tgen& gen)
        {
            vec X (d);
            for (auto& xk : X) xk = Z (gen);
            return value = X;
        }
};

class multi_normal: public random_variable<vec>
{
    private:
        mat T;
        normal Z;
    public:
        multi_normal (mat const& var_covar): T (chol (var_covar, "lower")), Z (0, 1) {};
        template<typename Tgen>
        vec operator () (Tgen& gen)
        {
            vec X (size (T.col (0)));
            for (auto& xk : X) xk = Z (gen);
            return value = T * X;
        }
};

#endif
