#pragma once

#include "../goldilocks/src/goldilocks_base_field.hpp"
#include "goldilocks_quadratic_ext.h"
#include <string>

//  改成点值存储
// class MultilinearPolynomial {
// public:
//     MultilinearPolynomial(size_t num_vars);
//     MultilinearPolynomial(const std::vector<Goldilocks2::Element>& evaluations);
//     // MultilinearPolynomial(){coeffs = {}; num_vars = 0;}
//     size_t get_num_vars() const{return num_vars;}
//     void set_term(std::string mask, const Goldilocks2::Element& c);
//     void set_term(std::string mask, const uint64_t& c);
//     Goldilocks2::Element get_term(uint64_t mask) const;
//     Goldilocks2::Element evaluate(const std::vector<Goldilocks2::Element>& point) const;
//     // MultilinearPolynomial operator+(const MultilinearPolynomial& other) const;
//     // MultilinearPolynomial operator*(const MultilinearPolynomial& other) const;
//     // static MultilinearPolynomial interpolate(const std::vector<Goldilocks2::Element>& values);
// private:
//     std::vector<Goldilocks2::Element> coeffs;
//     size_t num_vars;
// };


// store a multilinear polynomial in a vector of evaluations
class MultilinearPolynomial{
public:
    MultilinearPolynomial(size_t num_vars);
    MultilinearPolynomial(const std::vector<Goldilocks2::Element>& evaluations);
    size_t get_num_vars() const{return num_vars;}
    void set_value(const std::string& mask, const Goldilocks2::Element& c);
    void set_value(const std::string& mask, const uint64_t& c);
    Goldilocks2::Element eval_hypercube(uint64_t mask) const;
    Goldilocks2::Element evaluate(const std::vector<Goldilocks2::Element>& point) const;
    std::vector<Goldilocks2::Element> get_eval_table() const{return evaluations;}
    MultilinearPolynomial operator+(const MultilinearPolynomial& g) const;
    MultilinearPolynomial operator-(const MultilinearPolynomial& g) const;
private:
    // evaliations over the hypercube
    std::vector<Goldilocks2::Element> evaluations;
    size_t num_vars;
};
