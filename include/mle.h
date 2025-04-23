#pragma once

#include "../goldilocks/src/goldilocks_base_field.hpp"
#include "goldilocks_quadratic_ext.h"
#include <string>
class MultilinearPolynomial {
public:
    // using Goldilocks2::Element = Goldilocks2::Element;
    // using Goldilocks2::Element = Goldilocks2::Element;
    size_t num_vars;
    MultilinearPolynomial(size_t num_vars);

    void set_term(std::string mask, const Goldilocks2::Element& c);
    void set_term(std::string mask, const uint64_t& c);

    Goldilocks2::Element get_term(uint64_t mask) const;

    Goldilocks2::Element evaluate(const std::vector<Goldilocks2::Element>& point) const;

    // MultilinearPolynomial operator+(const MultilinearPolynomial& other) const;

    // MultilinearPolynomial operator*(const MultilinearPolynomial& other) const;

    // static MultilinearPolynomial interpolate(const std::vector<Goldilocks2::Element>& values);

private:
    std::vector<Goldilocks2::Element> coeffs;
    // std::vector<Goldilocks2::Element> keepTable;
};
