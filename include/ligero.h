#pragma once

#include "mle.h"
#include <cstdint>
#include <vector>

// class Ligero {
// public:
//     // commit a matrix, here w's evaluations instead of its coefficients are commited 
//     static inline std::vector<std::vector<Goldilocks2::Element>> commit(const MultilinearPolynomial& w);

// private:
//     double rho = 0.5;
// };
class ligeroProver{
public:
    ligeroProver(const MultilinearPolynomial& w, const double& rho);
    std::vector<std::vector<Goldilocks2::Element>> commit(const MultilinearPolynomial& w);
private:
    double rho;
    std::vector<Goldilocks2::Element> evals;
    size_t n;
    std::vector<std::vector<Goldilocks2::Element>> codewords;
};

class ligeroVerifier{
public:
    // std::vector<Goldilocks2::Element> commit(const MultilinearPolynomial& w);
    static inline bool verify(std::vector<Goldilocks2::Element> &Mhat, ligeroProver &prover);
    static inline Goldilocks2::Element open(ligeroProver &prover, const std::vector<Goldilocks2::Element> &r);
};

std::vector<Goldilocks2::Element> rsencode(const std::vector<Goldilocks2::Element> &data, const double& rho);