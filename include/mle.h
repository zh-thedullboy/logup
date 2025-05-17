#pragma once

#include "goldilocks_quadratic_ext.h"
#include <string>
#include <cstdint>


// store a multilinear polynomial in a vector of evaluations
class MultilinearPolynomial{
public:
    MultilinearPolynomial(size_t num_vars);
    MultilinearPolynomial(const std::vector<Goldilocks2::Element>& evaluations);
    MultilinearPolynomial(const std::vector<uint64_t>& val_table);
    size_t get_num_vars() const{return num_vars;}

    // only for debug use
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
