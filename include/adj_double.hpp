/* *********************************************************************************************** */
/*  COPYRIGHT (C) 2019 / 2020 SORBONNE UNIVERSITE - MSC PROBABILITY & FINANCE (EX - DEA EL KAROUI) */                                                                          
/*  ALL RIGHTS RESERVED. NO PART OF THIS PROJECT MAY BE REPRODUCED OR TRANSMITTED IN ANY FORM      */
/*  OR FOR ANY PURPOSE WITHOUT THE EXPRESS PERMISSION OF THE MANAGERS OF THE MASTER PROGRAM        */
/* *********************************************************************************************** */

#ifndef ADJDOUBLE_HPP
#define ADJDOUBLE_HPP

#include <vector>
#include <cmath>
#include <iostream>

#include "std_double.hpp"

struct Node
{
    uvec p; 
    vec w;

    Node () = default;
    Node (uvec const& p, vec const& w): p (p), w (w) {};
};

struct Tape
{
    std::vector<Node*> nodes;

    Tape () = default;
    Tape (Tape* t): nodes (t->nodes) {};

    unsigned size () const { return nodes.size (); };

    unsigned push0 ()
    {
        unsigned len (nodes.size ());
        nodes.push_back (new Node ());
        return len;
    }
    unsigned push (uvec const& p, vec const& w)
    {
        unsigned len (nodes.size ());
        nodes.push_back (new Node (p, w));
        return len;
    }

    template<typename adj_double>
    adj_double make0 (double value)
    {
        return adj_double (this, this->push0 (), value);
    }
    template<typename adj_double>
    vector<adj_double> make1 (vec const& values)
    {
        unsigned d (values.n_elem);
        vector<adj_double> u (d);
        for (unsigned i (0); i < d; ++i) u (i) = make0<adj_double> (values (i));
        return u;
    }
    template<typename adj_double>
    matrix<adj_double> make2 (mat const& values)
    {
        unsigned n (values.n_cols);
        matrix<adj_double> u (n);
        for (unsigned j (0); j < n; ++j) u (j) = make1<adj_double> (values.col (j));
        return u;
    }
    template<typename adj_double>
    matrix<adj_double> make_correlation_matrix (mat const& values)
    {
        unsigned n (values.n_cols);
        matrix<adj_double> u (n);
        for (unsigned j (0); j < n; ++j)
        {
            u (j) = vector<adj_double> (n);
            for (unsigned i (0); i < j; ++i) u (j) (i) = make0<adj_double> (values (i, j));
        } 
        return u;
    }
};

class Grad 
{
    public:
        Grad (vec const& derivatives): derivatives (derivatives) {};

        template<typename adj_double>
        double wrt (adj_double const& x) const { return derivatives (x.get_index ()); };

    private:
        vec derivatives;
};

class adj_double
{
    public:
        Tape* tape;
        unsigned index;
        double value;

    public:
        adj_double () = default;
        adj_double (Tape* tape, unsigned index, double value): tape (tape), index (index), value (value) {};

        void set_tape (Tape* t) { this->tape = t; };

        double get_value () const { return value; };
        unsigned get_index () const { return index; };

        adj_double operator-= (adj_double const& other)
        {
            value -= other.value;
            index = tape->push ({index, other.index}, {1, -1});
            return *this;
        }
        adj_double operator+= (adj_double const& other)
        {
            value += other.value;
            index = tape->push ({index, other.index}, {1, 1});
            return *this;
        }

        inline friend adj_double operator- (adj_double const& u, double c)
        {
            return adj_double (u.tape, u.tape->push ({u.index}, {1}), u.value - c);
        }
        inline friend adj_double operator+ (adj_double const& u, double c)
        {
            return adj_double (u.tape, u.tape->push ({u.index}, {1}), u.value + c);
        }
        inline friend adj_double operator* (adj_double const& u, adj_double const& v)
        {
            return adj_double (u.tape, u.tape->push ({u.index, v.index}, {v.value, u.value}), u.value * v.value);
        }
        inline friend adj_double operator* (adj_double const& u, double k)
        {
            return adj_double (u.tape, u.tape->push ({u.index}, {k}), k * u.value);
        }
        inline friend adj_double operator* (double k, adj_double const& u)
        {
            return u * k;
        }
        friend adj_double operator/ (adj_double const& u, adj_double const& v)
        {
            double inv_v = 1. / v.value;
            double u_inv_v = u.value * inv_v;
            double u_inv_sq_v = u_inv_v * inv_v;
            return adj_double (u.tape, u.tape->push ({u.index, v.index}, {inv_v, -u_inv_sq_v}), u_inv_v);
        }
        inline friend adj_double operator- (adj_double const& u)
        {
            return adj_double (u.tape, u.tape->push ({u.index}, {-1}), -u.value);
        }
        friend adj_double exp (adj_double const& u)
        {
            double exp_u (exp (u.value));
            return adj_double (u.tape, u.tape->push ({u.index}, {exp_u}), exp_u);
        }
        friend adj_double sqrt (adj_double const& u)
        {
            double sqrt_u (sqrt (u.value));
            return adj_double (u.tape, u.tape->push ({u.index}, {0.5 / sqrt_u}), sqrt_u);
        }
        friend adj_double max (adj_double const& u, double c)
        {
            double der (0), val (c);
            if (u.value > c) { val = u.value; der = 1; }
            return adj_double (u.tape, u.tape->push ({u.index}, {der}), val);
        }
        friend adj_double dot (vec const& w, vector<adj_double> const& u)
        {
            unsigned d (w.n_elem);
            uvec p (d);
            double sum (0);
            for (unsigned i (0); i < d; ++i) 
            {
                sum += w (i) * u (i).value;
                p (i) = u (i).index;
            }
            return adj_double (u (0).tape, u (0).tape->push (p, w), sum);
        }
        friend adj_double max (vector<adj_double> const& u)
        {
            unsigned d (u.n_elem);
            unsigned index_max (0);
            double max (-1e6), val (0);
            vec w (d, arma::fill::zeros);
            uvec p (d);
            for (unsigned i (0); i < d; ++i)
            {
                p (i) = u (i).index;
                val = u (i).value;
                if (val > max) { index_max = i; max = val; }
            } 
            w (index_max) = 1;
            return adj_double (u (0).tape, u (0).tape->push (p, w), max);
        }
        friend vector<adj_double> max (matrix<adj_double> const& U)
        {
            unsigned d (U.n_elem);
            vector<adj_double> u (d);
            for (unsigned i (0); i < d; ++i) u (i) = max (U (i));
            return u;
        }
        friend adj_double mean (vector<adj_double> const& u)
        {
            unsigned d (u.n_elem);
            double sum (0);
            vec w (d);
            uvec p (d);
            for (unsigned i (0); i < d; ++i)
            {
                w (i) = 1. / d;
                p (i) = u (i).index;
                sum += u (i).value;
            }
            return adj_double (u (0).tape, u (0).tape->push (p, w), sum / d);
        }
        friend void adj_cholesky (matrix<adj_double>& M, Tape* t)
        {
            unsigned n (M.n_elem);
            for (unsigned j = 0; j < n; ++j) 
            {
                auto l1 = t->make0<adj_double> (1.);
                for (unsigned k = 0; k < j; k++) l1 -= M (k) (j) * M (k) (j);
                M (j) (j) = sqrt (l1);
                for (unsigned i = j+1; i < n; ++i)
                {
                    auto l2 = M (i) (j); 
                    for (unsigned k = 0; k < j; k++) l2 -= M (k) (j) * M (k) (i);
                    M (j) (i) = l2 / M (j) (j);   
                }
            }
        }
        friend vector<adj_double> operator* (matrix<adj_double> const& M, vec const& x)
        {
            unsigned n (x.n_elem);
            vector<adj_double> y (n);
            for (unsigned i = 0; i < n; ++i)
            {
                auto sum = x (0) * M (0) (i);
                for (unsigned k (0); k < i; ++k) sum += x (k) * M (k) (i);
                y (i) = sum;
            }
            return y;
        }
        friend vector<adj_double> operator% (vec const& x, vector<adj_double> const& u)
        {
            unsigned n (x.n_elem);
            vector<adj_double> v (n);
            for (unsigned i = 0; i < n; ++i) v (i) = x (i) * u (i);
            return v;
        }
        friend vector<adj_double> exp (vector<adj_double> const& u)
        {
            unsigned n (u.n_elem);
            vector<adj_double> v (n);
            for (unsigned i = 0; i < n; ++i) v (i) = exp (u (i));
            return v;
        }
        friend vector<adj_double> operator+ (vector<adj_double> const& u, vec const& x)
        {
            unsigned n (x.n_elem);
            vector<adj_double> v (n);
            for (unsigned i = 0; i < n; ++i) v (i) = u (i) + x (i);
            return v;
        }
        friend vector<adj_double> operator* (double k, vector<adj_double> const& x)
        {
            unsigned n (x.n_elem);
            vector<adj_double> y (n);
            for (unsigned i = 0; i < n; ++i) y (i) = k * x (i);
            return y;
        }

        Grad grad (double seed = 1) const
        {
            vec derivatives (tape->size (), arma::fill::zeros);
            auto nodes = tape->nodes;
            derivatives (index) = seed;
            for (unsigned i (index); i > 0; i--)
            {
                auto node = nodes [i];
                auto derivative = derivatives (i);
                auto indexes = node->p;
                auto weights = node->w;
                auto it_indexes = indexes.begin ();
                auto it_weights = weights.begin ();
                if (!indexes.empty ())
                    while (it_indexes != indexes.end () && it_weights != weights.end ())
                        derivatives (*it_indexes++) += (*it_weights++) * derivative;
            }
            return Grad (derivatives);
        }
};

#endif

// Tree-based implementation

/*

template<typename T>
struct Child
{
    using Node = T;
    private:
        double child_arc_derivative;
        Node child_node;
    public:
        Child () = default;
        Child (double child_arc_derivative, Node const& child_node): 
            child_arc_derivative (child_arc_derivative), child_node (child_node) {}
        double get_child_arc_derivative () { return this->child_arc_derivative; }
        Node get_child_node () { return this->child_node; }
};

struct Node
{
    private:
        double value;
        std::vector <Child <Node>> *children;
        double adjoint;
    public:
        Node () = default;
        ~Node () = default;
        Node (double value): value (value), children (new std::vector<Child<Node>>), adjoint (1e-50) {}
        friend std::ostream& operator<< (std::ostream& os, Node const& x) 
            { return os << x.value << std::endl; }
        void clear () { this->children->clear (); }
        double get_adjoint (double lambda = 1)
        {
            if (this->adjoint == 1e-50)
            {
                if (this->children->empty ()) return lambda;
                else 
                {   
                    double sum (0);
                    for (auto& e : *(this->children)) 
                        sum += e.get_child_arc_derivative () * e.get_child_node ().get_adjoint (lambda);
                    return this->adjoint = sum;
                }
            } else return this->adjoint;
        }
        friend Node operator+ (Node const& x, Node const& y)
        {
            Node z (x.value + y.value);
            x.children->push_back (Child<Node> (1, z));
            y.children->push_back (Child<Node> (1, z));
            return z;
        }
        friend Node operator- (Node const& x, Node const& y)
        {
            Node z (x.value - y.value);
            x.children->push_back (Child<Node> (1, z));
            y.children->push_back (Child<Node> (-1, z));
            return z;
        }
        friend Node operator* (Node const& x, Node const& y)
        {
            Node z (x.value * y.value);
            x.children->push_back (Child<Node> (y.value, z));
            y.children->push_back (Child<Node> (x.value, z));
            return z;
        }
        friend Node operator/ (Node const& x, Node const& y)
        {
            double inv_y_value (1. / y.value);
            Node z (x.value * inv_y_value);
            x.children->push_back (Child<Node> (inv_y_value, z));
            y.children->push_back (Child<Node> (-z.value * inv_y_value, z));
            return z;
        }
        friend Node cos (Node const& x)
        {
            Node z (cos (x.value));
            x.children->push_back (Child<Node> (-sin (x.value), z));
            return z;
        }
        friend Node max (Node const& x, Node const& y)
        {
            double max (x.value), dx (1), dy (0);
            if (x.value < y.value) { max = y.value; dx = 0; dy = 1; }
            Node z (max);
            x.children->push_back (Child<Node> (dx, z));
            y.children->push_back (Child<Node> (dy, z));
            return z;
        }
        friend Node exp (Node const& x)
        {
            double exp_x (exp (x.value));
            Node z (exp_x);
            x.children->push_back (Child<Node> (exp_x, z));
            return z;
        }
        friend Node log (Node const& x)
        {
            Node z (log (x.value));
            x.children->push_back (Child<Node> (1. / x.value, z));
            return z;
        }
        friend Node pow (Node const& x, unsigned n)
        {
            double pow_1 (pow (x.value, n - 1));
            Node z (x.value * pow_1);
            x.children->push_back (Child<Node> (n * pow_1, z));
            return z;
        }
};

*/