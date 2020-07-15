/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#ifndef TGTDOUBLE_HPP
#define TGTDOUBLE_HPP

#include <cmath>
#include <iostream>

#include "toolBox.hpp"

class tgt_double
{
    using value_type = double;
    using derivative_type = vec;

    private:

        value_type f; 
        derivative_type df;

    public:

        tgt_double () = default;
        tgt_double (value_type const& f, derivative_type const& df): f (f), df (df) {};

        value_type get_value () const { return f; };
        derivative_type get_derivative () const { return df; };
        
        tgt_double operator+= (tgt_double const& v)
        {
            f += v.f;
            df += v.df;
            return *this;
        }
        inline friend tgt_double operator+ (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f + v.f, u.df + v.df);
        }
        inline friend tgt_double operator/ (tgt_double const& u, double k)
        {
            return tgt_double (u.f / k, u.df / k);
        }
        inline friend tgt_double operator- (tgt_double const& u, double c)
        {
            return tgt_double (u.f - c, u.df);
        }
        inline friend tgt_double operator- (tgt_double const& u)
        {
            return tgt_double (-u.f, -u.df);
        }
        inline friend tgt_double operator* (double k, tgt_double const& u)  
        {
            return tgt_double (k * u.f, k * u.df);
        }
        inline friend tgt_double operator* (tgt_double const& u, double k)  
        {
            return k * u;
        }
        inline friend tgt_double operator* (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f * v.f, v.f * u.df + u.f * v.df);
        }
        inline friend tgt_double operator/ (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f / v.f, (v.f * u.df - u.f * v.df) / v.f / v.f);
        }
        inline friend tgt_double operator/ (double k, tgt_double const& v)
        {
            return tgt_double (k / v.f, -k * v.df / v.f / v.f);
        }
        friend tgt_double exp (tgt_double const& u)
        {
            double exp_u = exp (u.f);
            return tgt_double (exp_u, exp_u * u.df);
        }
        inline friend tgt_double log (tgt_double const& u)
        {
            return tgt_double (log (u.f), u.df / u.f);
        }
        inline friend tgt_double cos (tgt_double const& u)
        {
            return tgt_double (cos (u.f), -sin (u.f) * u.df);
        }
        inline friend tgt_double sin (tgt_double const& u)
        {
            return tgt_double (sin (u.f), cos (u.f) * u.df);
        }
        inline friend tgt_double tan (tgt_double const& u)
        {
            return tgt_double (tan (u.f), u.df / (1 + u.f * u.f));
        }
        inline friend tgt_double max (tgt_double const& u, value_type const& c)
        {
            if (u.f > c) return tgt_double (u.f, u.df); else return tgt_double (c, arma::zeros (size (u.df)));
        }
        friend tgt_double dot (vec const& w, vector<tgt_double> const& u)
        {
            unsigned d (w.n_elem);
            auto sum = w (0) * u (0); 
            for (unsigned i (1); i < d; ++i) sum += w (i) * u (i);
            return sum;
        }
        friend tgt_double max (vector<tgt_double> const& u)
        {
            unsigned n (u.n_elem);
            unsigned index_max (0);
            double max (-1e6), val (0); 
            vec der (size (u (0).df));
            for (unsigned i (0); i < n; ++i)
            {
                val = u (i).f;
                if (val > max) { max = val; index_max = i;}
            } 
            return u (index_max);
        }
        friend vector<tgt_double> max (matrix<tgt_double> const& U)
        {
            unsigned n (U.n_elem);
            vector<tgt_double> u (n);
            for (unsigned i (0); i < n; ++i) u (i) = max (U (i));
            return u;
        }
        friend tgt_double mean (vector<tgt_double> const& u)
        {
            unsigned n (u.n_elem);
            auto sum (u (0));
            for (unsigned i (1); i < n; ++i) sum += u (i);
            return sum / n;
        }
};

template<typename tgt_double>
tgt_double make0 (double f, vec const& df)
{
    return tgt_double (f, df);
}

template<typename tgt_double>
vector<tgt_double> make1 (vec const& f, mat const& df)
{
    unsigned n (f.n_elem);
    vector<tgt_double> res (n);
    for (unsigned i (0); i < n; ++i) res (i) = make0<tgt_double> (f (i), df.col (i));
    return res;
}

template<typename tgt_double>
matrix<tgt_double> make2 (mat const& f, cube const& df)
{
    unsigned m (f.n_cols);
    matrix<tgt_double> res (m);
    for (unsigned i (0); i < m; ++i) res (i) = make1<tgt_double> (f.col (i), df.slice (i));
    return res;
}

#endif

/*

template<unsigned d>
using vec = arma::vec::fixed<d>;

template<unsigned n, unsigned d>
using mat = arma::mat::fixed<n, d>;

template<unsigned n>
class tgt_double
{
    using value_type = double;
    using derivative_type = vec<n>;
    private:
        value_type f; 
        derivative_type df;
    public:
        tgt_double () = default;
        ~tgt_double () = default;
        tgt_double (value_type const& f, derivative_type const& df): f (f), df (df) {};
        value_type get_value () const { return f; };
        derivative_type get_derivative () const { return df; };
        friend std::ostream& operator<< (std::ostream& os, tgt_double const& x)
        {
            return os << 
                "Value: " << x.f << std::endl <<
                "Derivative: " << std::endl << x.df << std::endl;
        }       
        inline friend tgt_double operator+ (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f + v.f, u.df + v.df);
        }
        inline friend tgt_double operator+ (tgt_double const& u, double const& c) 
        { 
            return tgt_double (u.f + c, u.df); 
        }
        tgt_double operator+= (tgt_double const& u) 
        { 
            f += u.f;
            df += u.df;
            return *this;
        }
        inline friend tgt_double sqrt (tgt_double const& u) 
        { 
            value_type sqrt_u (sqrt (u.f));
            return tgt_double (sqrt_u, u.df / sqrt_u); 
        }
        inline friend tgt_double operator% (tgt_double const& u, tgt_double const& v) 
        { 
            return tgt_double (u.f * v.f, u.df % v.df); 
        }
        inline friend tgt_double operator- (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f - v.f, u.df - v.df);
        }
        inline friend tgt_double operator- (tgt_double const& u, double const& c)
        {
            return tgt_double (u.f - c, u.df);
        }
        inline friend tgt_double operator* (double const& k, tgt_double const& u)  
        {
            return tgt_double (k * u.f, k * u.df);
        }
        inline friend tgt_double operator* (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f * v.f, v.f * u.df + u.f * v.df);
        }
        inline friend tgt_double operator/ (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f / v.f, (v.f * u.df - u.f * v.df) / v.f / v.f);
        }
        inline friend tgt_double operator/ (double const& k, tgt_double const& v)
        {
            return tgt_double (k / v.f, -k * v.df / v.f / v.f);
        }
        inline friend tgt_double operator/ (tgt_double const& u, double const& k)
        {
            return tgt_double (u.f / k, u.df / k);
        }
        inline friend tgt_double sin (tgt_double const& u) 
        {
            return tgt_double (sin (u.f), cos (u.f) * u.df);
        }
        inline friend tgt_double cos (tgt_double const& u)
        {
            return tgt_double (cos (u.f), -sin (u.f) * u.df);
        }
        inline friend tgt_double tan (tgt_double const& u)
        {
            value_type tan_u = tan (u.f);
            return tgt_double (tan_u, (1 + tan_u * tan_u) * u.df);
        }
        inline friend tgt_double pow (tgt_double const& u, unsigned p)
        {
            double pow_u = pow (u.f, p - 1);
            return tgt_double (pow_u * u.f, pow_u * u.df);
        }
        inline friend tgt_double exp (tgt_double const& u)
        {
            double exp_u = exp (u.f);
            return tgt_double (exp_u, exp_u * u.df);
        }
        inline friend tgt_double log (tgt_double const& u)
        {
            return tgt_double (log (u.f), u.df / u.f);
        }
        inline friend tgt_double max (tgt_double const& u, double const& c)
        {
            if (u.f > c) return tgt_double (u.f, u.df); else return tgt_double (c, arma::zeros (n));
        }
};

template<unsigned d, unsigned n> 
class bigtgt_double
{
    friend class tgt_double<n>;
    using value_type = vec<d>;
    using derivative_type = mat<n, d>;
    private:
        value_type f;
        derivative_type df; 
    public:
        bigtgt_double () = default;
        ~bigtgt_double () = default;
        bigtgt_double (value_type const& f, derivative_type const& df): f (f), df (df) {};
        value_type get_value () const { return f; };
        derivative_type get_derivative () const { return df; };
        friend std::ostream& operator<< (std::ostream& os, bigtgt_double const& x)
        {
            return os << 
                "Value: " << std::endl << x.f <<
                "Derivative: " << std::endl << x.df;
        }       
        inline friend bigtgt_double operator+ (bigtgt_double const& u, bigtgt_double const& v)
        {
            return bigtgt_double (u.f + v.f, u.df + v.df);
        }
        inline friend bigtgt_double operator+ (bigtgt_double const& u, value_type const& c) 
        { 
            return bigtgt_double (u.f + c, u.df); 
        }
        inline friend bigtgt_double operator- (bigtgt_double const& u, bigtgt_double const& v)
        {
            return bigtgt_double (u.f - v.f, u.df - v.df);
        }
        inline friend bigtgt_double operator- (bigtgt_double const& u, value_type const& c)
        {
            return bigtgt_double (u.f - c, u.df);
        }
        inline friend bigtgt_double operator* (double const& k, bigtgt_double const& u)  
        {
            return bigtgt_double (k * u.f, k * u.df);
        }
        inline friend bigtgt_double operator/ (bigtgt_double const& u, double const& k)
        {
            return bigtgt_double (u.f / k, u.df / k);
        }
        inline friend tgt_double<n> dot (bigtgt_double const& u, bigtgt_double const& v)
        {
            return tgt_double<n> (dot (u.f, v.f), v.df * u.f + u.df * v.f);
        }
        inline friend tgt_double<n> dot (bigtgt_double const& u, value_type const& w)
        {
            return tgt_double<n> (dot (u.f, w), u.df * w);
        }
        inline friend bigtgt_double operator/ (bigtgt_double const& u, bigtgt_double const& v)
        {
            auto u_ = u.df;
            auto v_ = v.df;
            auto sq_v = v.f % v.f;
            u_.each_row () %= (v.f / sq_v).as_row ();
            v_.each_row () %= (u.f / sq_v).as_row ();
            return bigtgt_double (u.f / v.f, u_ - v_);
        }
        inline friend bigtgt_double operator/ (bigtgt_double const& u, value_type const& k)
        {
            auto u_ = u.df;
            u_.each_row () /= k.as_row ();
            return bigtgt_double (u.f / k, u_);
        }
        inline friend bigtgt_double operator/ (value_type const& k, bigtgt_double const& v)
        {
            auto v_ = v.df;
            auto sq_v = v.f % v.f;
            v_.each_row () %= -(k / sq_v).as_row ();
            return bigtgt_double (k / v.f, v_);
        }
        inline friend bigtgt_double sin (bigtgt_double const& u) 
        {
            auto u_ = u.df;
            u_.each_row () %= cos (u.f).as_row ();
            return bigtgt_double (sin (u.f), u_);
        }
        inline friend bigtgt_double cos (bigtgt_double const& u)
        {
            auto u_ = u.df;
            u_.each_row () %= -sin (u.f).as_row ();
            return bigtgt_double (cos (u.f), u_);
        }
        inline friend bigtgt_double exp (bigtgt_double const& u)
        {
            auto u_ = u.df;
            value_type exp_u = exp (u.f);
            u_.each_row () %= exp_u.as_row ();
            return bigtgt_double (exp_u, u_);
        }
        inline friend bigtgt_double log (bigtgt_double const& u)
        {
            auto u_ = u.df;
            u_.each_row () /= (u.f).as_row ();
            return bigtgt_double (log (u.f), u_);
        }
        inline friend bigtgt_double operator% (bigtgt_double const& u, bigtgt_double const& v)
        {
            auto u_ = u.df;
            auto v_ = v.df;
            u_.each_row () %= (v.f).as_row ();
            v_.each_row () %= (u.f).as_row ();
            return bigtgt_double (u.f % v.f, u_ + v_);
        }
        inline friend bigtgt_double operator% (value_type const& k, bigtgt_double const& u)
        {
            auto u_ = u.df;
            u_.each_row () %= k.as_row ();
            return bigtgt_double (k % u.f, u_);
        }
};

*/

/*

using vec = arma::vec;
using mat = arma::mat;

template<typename T1, typename T2>
class tgt_double
{
    using value_type = T1;
    using derivative_type = T2;
    public:
        value_type f; 
        derivative_type df;
        tgt_double () = default;
        tgt_double (value_type const& f, derivative_type const& df): f (f), df (df) {};
        ~tgt_double () = default;
        value_type get_value () const { return f; };
        derivative_type get_derivative () const { return df; };
        friend std::ostream& operator<< (std::ostream& os, tgt_double const& x)
        {
            return os << "Value: " << x.f << std::endl << "Derivative: " << x.df << std::endl;
        }
        inline friend tgt_double operator- (tgt_double const& u, value_type const& c)
        {
            return tgt_double (u.f - c, u.df);
        }
        inline friend tgt_double operator- (tgt_double const& u)
        {
            return tgt_double (-u.f, -u.df);
        }
        inline friend tgt_double operator* (double k, tgt_double const& u)  
        {
            return tgt_double (k * u.f, k * u.df);
        }
        inline friend tgt_double operator* (tgt_double const& u, double k)  
        {
            return k * u;
        }     
};

template<typename T1, typename T2>
tgt_double<T1,T2> make_tgt_double (T1 const& value, T2 const& derivative)
{
    return tgt_double<T1,T2> (value, derivative);
}

tgt_double<double, vec> operator* (tgt_double<double, vec> const& u, tgt_double<double, vec> const& v)
{
    return tgt_double<double, vec> (u.f * v.f, v.f * u.df + u.f * v.df);
}
tgt_double<double, vec> exp (tgt_double<double, vec> const& u)
{
    double exp_u = exp (u.f);
    return tgt_double<double, vec> (exp_u, exp_u * u.df);
}
tgt_double<double, vec> max (tgt_double<double, vec> const& u, double c)
{
    if (u.f > c) return tgt_double<double, vec> (u.f, u.df); else return tgt_double<double, vec> (c, arma::zeros (size (u.df)));
}
tgt_double<double, vec> dot (vec const& w, tgt_double<vec, mat> const& u)
{
    return tgt_double<double, vec> (dot (w, u.f), u.df * w);
}
    [...]



*/

/*


#include "array.hpp"

class tgt_double
{
    using value_type = double;
    using derivative_type = arma::vec;
    private:
        value_type f; 
        derivative_type df;
    public:
        tgt_double () = default;
        tgt_double (value_type const& f, derivative_type const& df): f (f), df (df) {};
        tgt_double (value_type const& f, unsigned n, unsigned i): f (f), df (n, 0) { df (i) = 1; };
        tgt_double (unsigned d): f (0), df (d) { df.fill (0); };
        ~tgt_double () = default;
        value_type get_value () const { return f; };
        derivative_type get_derivative () const { return df; };
        friend std::ostream& operator<< (std::ostream& os, tgt_double const& x)
        {
            return os << "Value: " << x.f << std::endl << "Derivative: " << x.df << std::endl;
        }
        tgt_double operator+= (tgt_double const& v)
        {
            f += v.f;
            df += v.df;
            return *this;
        }
        inline friend tgt_double operator+ (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f + v.f, u.df + v.df);
        }
        inline friend tgt_double operator+ (tgt_double const& u, double c) 
        { 
            return tgt_double (u.f + c, u.df); 
        }
        friend tgt_double sqrt (tgt_double const& u) 
        { 
            double sqrt_u (sqrt (u.f));
            return tgt_double (sqrt_u, 0.5 / sqrt_u * u.df); 
        }
        inline friend tgt_double operator- (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f - v.f, u.df - v.df);
        }
        inline friend tgt_double operator- (tgt_double const& u, double c)
        {
            return tgt_double (u.f - c, u.df);
        }
        inline friend tgt_double operator- (tgt_double const& u)
        {
            return tgt_double (-u.f, -u.df);
        }
        inline friend tgt_double operator* (double k, tgt_double const& u)  
        {
            return tgt_double (k * u.f, k * u.df);
        }
        inline friend tgt_double operator* (tgt_double const& u, double k)  
        {
            return k * u;
        }
        inline friend tgt_double operator* (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f * v.f, v.f * u.df + u.f * v.df);
        }
        inline friend tgt_double operator/ (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f / v.f, (v.f * u.df - u.f * v.df) / v.f / v.f);
        }
        inline friend tgt_double operator/ (double k, tgt_double const& v)
        {
            return tgt_double (k / v.f, -k * v.df / v.f / v.f);
        }
        inline friend tgt_double operator/ (tgt_double const& u, double k)
        {
            return tgt_double (u.f / k, u.df / k);
        }
        inline friend tgt_double sin (tgt_double const& u) 
        {
            return tgt_double (sin (u.f), cos (u.f) * u.df);
        }
        inline friend tgt_double cos (tgt_double const& u)
        {
            return tgt_double (cos (u.f), -sin (u.f) * u.df);
        }
        friend tgt_double tan (tgt_double const& u)
        {
            double tan_u = tan (u.f);
            return tgt_double (tan_u, (1 + tan_u * tan_u) * u.df);
        }
        friend tgt_double pow (tgt_double const& u, unsigned n)
        {
            double pow_u = pow (u.f, n - 1);
            return tgt_double (pow_u * u.f, n * pow_u * u.df);
        }
        friend tgt_double exp (tgt_double const& u)
        {
            double exp_u = exp (u.f);
            return tgt_double (exp_u, exp_u * u.df);
        }
        inline friend tgt_double log (tgt_double const& u)
        {
            return tgt_double (log (u.f), u.df / u.f);
        }
        inline friend tgt_double max (tgt_double const& u, double const& c)
        {
            if (u.f > c) return tgt_double (u.f, u.df); else return tgt_double (c, arma::zeros (arma::size (u.df)));
        }
        friend tgt_double dot (arma::vec const& w, array<tgt_double> const& u)
        {
            //auto n (arma::size (u(0).df));
            //auto d (arma::size (w));
            auto n (2);
            auto d (2);
            tgt_double sum (n); 
            for (unsigned i (0); i < d; ++i) sum += w (i) * u (i);
            return sum;
        }
};

template<unsigned d>
using vec = arma::vec::fixed<d>;

template<unsigned n, unsigned d>
using mat = arma::mat::fixed<n, d>;

template<unsigned n>
class tgt_double
{
    using value_type = double;
    using derivative_type = vec<n>;
    private:
        value_type f; 
        derivative_type df;
    public:
        tgt_double () = default;
        ~tgt_double () = default;
        tgt_double (value_type const& f, derivative_type const& df): f (f), df (df) {};
        value_type get_value () const { return f; };
        derivative_type get_derivative () const { return df; };
        friend std::ostream& operator<< (std::ostream& os, tgt_double const& x)
        {
            return os << 
                "Value: " << x.f << std::endl <<
                "Derivative: " << std::endl << x.df << std::endl;
        }       
        inline friend tgt_double operator+ (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f + v.f, u.df + v.df);
        }
        inline friend tgt_double operator+ (tgt_double const& u, double const& c) 
        { 
            return tgt_double (u.f + c, u.df); 
        }
        tgt_double operator+= (tgt_double const& u) 
        { 
            f += u.f;
            df += u.df;
            return *this;
        }
        inline friend tgt_double sqrt (tgt_double const& u) 
        { 
            value_type sqrt_u (sqrt (u.f));
            return tgt_double (sqrt_u, u.df / sqrt_u); 
        }
        inline friend tgt_double operator% (tgt_double const& u, tgt_double const& v) 
        { 
            return tgt_double (u.f * v.f, u.df % v.df); 
        }
        inline friend tgt_double operator- (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f - v.f, u.df - v.df);
        }
        inline friend tgt_double operator- (tgt_double const& u, double const& c)
        {
            return tgt_double (u.f - c, u.df);
        }
        inline friend tgt_double operator* (double const& k, tgt_double const& u)  
        {
            return tgt_double (k * u.f, k * u.df);
        }
        inline friend tgt_double operator* (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f * v.f, v.f * u.df + u.f * v.df);
        }
        inline friend tgt_double operator/ (tgt_double const& u, tgt_double const& v)
        {
            return tgt_double (u.f / v.f, (v.f * u.df - u.f * v.df) / v.f / v.f);
        }
        inline friend tgt_double operator/ (double const& k, tgt_double const& v)
        {
            return tgt_double (k / v.f, -k * v.df / v.f / v.f);
        }
        inline friend tgt_double operator/ (tgt_double const& u, double const& k)
        {
            return tgt_double (u.f / k, u.df / k);
        }
        inline friend tgt_double sin (tgt_double const& u) 
        {
            return tgt_double (sin (u.f), cos (u.f) * u.df);
        }
        inline friend tgt_double cos (tgt_double const& u)
        {
            return tgt_double (cos (u.f), -sin (u.f) * u.df);
        }
        inline friend tgt_double tan (tgt_double const& u)
        {
            value_type tan_u = tan (u.f);
            return tgt_double (tan_u, (1 + tan_u * tan_u) * u.df);
        }
        inline friend tgt_double pow (tgt_double const& u, unsigned p)
        {
            double pow_u = pow (u.f, p - 1);
            return tgt_double (pow_u * u.f, pow_u * u.df);
        }
        inline friend tgt_double exp (tgt_double const& u)
        {
            double exp_u = exp (u.f);
            return tgt_double (exp_u, exp_u * u.df);
        }
        inline friend tgt_double log (tgt_double const& u)
        {
            return tgt_double (log (u.f), u.df / u.f);
        }
        inline friend tgt_double max (tgt_double const& u, double const& c)
        {
            if (u.f > c) return tgt_double (u.f, u.df); else return tgt_double (c, arma::zeros (n));
        }
};

template<unsigned d, unsigned n> 
class bigtgt_double
{
    friend class tgt_double<n>;
    using value_type = vec<d>;
    using derivative_type = mat<n, d>;
    private:
        value_type f;
        derivative_type df; 
    public:
        bigtgt_double () = default;
        ~bigtgt_double () = default;
        bigtgt_double (value_type const& f, derivative_type const& df): f (f), df (df) {};
        value_type get_value () const { return f; };
        derivative_type get_derivative () const { return df; };
        friend std::ostream& operator<< (std::ostream& os, bigtgt_double const& x)
        {
            return os << 
                "Value: " << std::endl << x.f <<
                "Derivative: " << std::endl << x.df;
        }       
        inline friend bigtgt_double operator+ (bigtgt_double const& u, bigtgt_double const& v)
        {
            return bigtgt_double (u.f + v.f, u.df + v.df);
        }
        inline friend bigtgt_double operator+ (bigtgt_double const& u, value_type const& c) 
        { 
            return bigtgt_double (u.f + c, u.df); 
        }
        inline friend bigtgt_double operator- (bigtgt_double const& u, bigtgt_double const& v)
        {
            return bigtgt_double (u.f - v.f, u.df - v.df);
        }
        inline friend bigtgt_double operator- (bigtgt_double const& u, value_type const& c)
        {
            return bigtgt_double (u.f - c, u.df);
        }
        inline friend bigtgt_double operator* (double const& k, bigtgt_double const& u)  
        {
            return bigtgt_double (k * u.f, k * u.df);
        }
        inline friend bigtgt_double operator/ (bigtgt_double const& u, double const& k)
        {
            return bigtgt_double (u.f / k, u.df / k);
        }
        inline friend tgt_double<n> dot (bigtgt_double const& u, bigtgt_double const& v)
        {
            return tgt_double<n> (dot (u.f, v.f), v.df * u.f + u.df * v.f);
        }
        inline friend tgt_double<n> dot (bigtgt_double const& u, value_type const& w)
        {
            return tgt_double<n> (dot (u.f, w), u.df * w);
        }
        inline friend bigtgt_double operator/ (bigtgt_double const& u, bigtgt_double const& v)
        {
            auto u_ = u.df;
            auto v_ = v.df;
            auto sq_v = v.f % v.f;
            u_.each_row () %= (v.f / sq_v).as_row ();
            v_.each_row () %= (u.f / sq_v).as_row ();
            return bigtgt_double (u.f / v.f, u_ - v_);
        }
        inline friend bigtgt_double operator/ (bigtgt_double const& u, value_type const& k)
        {
            auto u_ = u.df;
            u_.each_row () /= k.as_row ();
            return bigtgt_double (u.f / k, u_);
        }
        inline friend bigtgt_double operator/ (value_type const& k, bigtgt_double const& v)
        {
            auto v_ = v.df;
            auto sq_v = v.f % v.f;
            v_.each_row () %= -(k / sq_v).as_row ();
            return bigtgt_double (k / v.f, v_);
        }
        inline friend bigtgt_double sin (bigtgt_double const& u) 
        {
            auto u_ = u.df;
            u_.each_row () %= cos (u.f).as_row ();
            return bigtgt_double (sin (u.f), u_);
        }
        inline friend bigtgt_double cos (bigtgt_double const& u)
        {
            auto u_ = u.df;
            u_.each_row () %= -sin (u.f).as_row ();
            return bigtgt_double (cos (u.f), u_);
        }
        inline friend bigtgt_double exp (bigtgt_double const& u)
        {
            auto u_ = u.df;
            value_type exp_u = exp (u.f);
            u_.each_row () %= exp_u.as_row ();
            return bigtgt_double (exp_u, u_);
        }
        inline friend bigtgt_double log (bigtgt_double const& u)
        {
            auto u_ = u.df;
            u_.each_row () /= (u.f).as_row ();
            return bigtgt_double (log (u.f), u_);
        }
        inline friend bigtgt_double operator% (bigtgt_double const& u, bigtgt_double const& v)
        {
            auto u_ = u.df;
            auto v_ = v.df;
            u_.each_row () %= (v.f).as_row ();
            v_.each_row () %= (u.f).as_row ();
            return bigtgt_double (u.f % v.f, u_ + v_);
        }
        inline friend bigtgt_double operator% (value_type const& k, bigtgt_double const& u)
        {
            auto u_ = u.df;
            u_.each_row () %= k.as_row ();
            return bigtgt_double (k % u.f, u_);
        }
};
*/