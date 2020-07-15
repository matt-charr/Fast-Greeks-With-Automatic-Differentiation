/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <cmath>

#include "armadillo.hpp"

template <typename T>
class mean_var
{
    using result_type = T;
    public:
        mean_var () = default;
        mean_var (unsigned n, T const& sumX, T const& sumXX): 
            sampleSize (n), sumX (sumX), sumXX (sumXX) {};
        result_type mean () const
        {
            return sumX / (double) sampleSize;
        };
        result_type var () const
        {
            result_type tmp (mean ());
            return sumXX / sampleSize - tmp % tmp;
        };
        result_type icSize () const
        {
            return 1.33 * sqrt ((var ()) / sampleSize);
        };
        mean_var& operator+= (mean_var const& mv)
        {
            sampleSize += mv.sampleSize;
            sumX += mv.sumX;
            sumXX += mv.sumXX;
            return *this;
        };
        friend mean_var operator+ (mean_var const& mv1, mean_var const& mv2)
        {
            mean_var mv (mv1.sampleSize + mv2.sampleSize, mv1.sumX + mv2.sumX, mv1.sumXX + mv2.sumXX);
            return mv; 
        };
        friend std::ostream& operator<< (std::ostream& o, mean_var const& mv)
        {
            return o << std::endl << "@Mean" << std::endl << std::endl << mv.mean () << std::endl << "#icSize" << std::endl << std::endl << mv.icSize () << std::endl;
        }
    private :
        unsigned sampleSize;
        result_type sumX, sumXX;
};

template <typename T, typename Tdistribution, typename Tgen>
mean_var<T> monteCarlo (Tdistribution& X, unsigned n, Tgen& gen)
{
    auto x = X (gen);
    auto sumX = x;
    auto sumXX = x; sumXX %= x;
    for(unsigned i (1); i < n; ++i)
    {
        x = X (gen);
        sumX += x;
        sumXX += x % x;
    }
    return mean_var<T> (n, sumX, sumXX);
}

template<typename T, typename Tdistribution, typename Tgen>
float timeMonteCarlo (Tdistribution& X, unsigned n, Tgen& gen)
{
    clock_t t1, t2;
    t1 = clock ();
    monteCarlo<T, Tdistribution, Tgen> (X, n, gen);
    t2 = clock ();
    return (float) (t2-t1) / CLOCKS_PER_SEC;
}

template <typename Outer, typename Inner>
class composed
{
    private:
        Outer f;
        Inner g;
    public: 
        composed (Outer f, Inner g) : f(f), g(g) {}
        template<typename T>
        auto operator() (T& x) { return f (g (x)); };
};

template <typename Outer, typename Inner>
inline composed<Outer, Inner> compose (Outer f, Inner g)
{
    return composed<Outer, Inner> (f, g);
}

#endif