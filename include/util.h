#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include "goldilocks_quadratic_ext.h"
#include "mle.h"

uint64_t convert_mask_to_u64(const std::string& mask, const size_t &nvar);

int find_ceiling_log2(const uint64_t& n);

// evaluate \tilde{eq}(r, x) = \prod_{i=0}^{n-1} (1 - r_i x_i) in O(2^l) linear time
MultilinearPolynomial eq(const size_t& num_var, const std::vector<Goldilocks2::Element>& r);

// calculate the inverse of all elements in arr with calculating only one inverse
void batch_inverse(std::vector<Goldilocks2::Element>& inv, const std::vector<Goldilocks2::Element>& arr);

std::vector<Goldilocks2::Element> random_vec(size_t n);

void print_table(const std::vector<Goldilocks2::Element>& table);

void print_table(const std::vector<size_t>& table);

void print_hash(const std::array<uint8_t, 32>& hash);
void print_bytes(const std::array<uint8_t, 16>& hash);

// evaluate a polynomial with its coefficients known as coefs at point x with Horners method
Goldilocks2::Element Horner(const std::vector<Goldilocks2::Element> &coefs, const Goldilocks2::Element& x);

// MultilinearPolynomial buildpolynomial(const Goldilocks::Element& gamma, const Goldilocks::Element& lambda, std::vector<Goldilocks2::Element>& p1, std::vector<Goldilocks2::Element>& p2);

namespace std{
    template <>
    struct hash<Goldilocks2::Element> {
        size_t operator()(const Goldilocks2::Element& arr) const {
            size_t h1 = std::hash<uint64_t>{}(Goldilocks::toU64(arr[0]));
            size_t h2 = std::hash<uint64_t>{}(Goldilocks::toU64(arr[1]));
            return h1 ^ (h2 << 1);
        }
    };
}

