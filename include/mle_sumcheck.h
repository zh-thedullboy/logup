#pragma once

#include "mle.h"
#include "goldilocks_quadratic_ext.h"
#include <vector>
#include <array>

class Prover{
public:
    Prover(const MultilinearPolynomial& g);
    void initialize();
    std::array<Goldilocks2::Element, 2> send_message(const size_t& round,const std::vector<Goldilocks2::Element>& rands);
    Goldilocks2::Element get_sum() const { return sum; }
    size_t get_rounds() const { return nrnd; }
    private:
    MultilinearPolynomial g;
    std::vector<Goldilocks2::Element> keepTable;
    Goldilocks2::Element sum;
    size_t nrnd;
};

class Verifier{
public:
    // should be replaced with a pcs
    typedef MultilinearPolynomial Oracle;
    static bool execute_sumcheck(Prover& pr, const Oracle& oracle);
private:
    static Goldilocks2::Element challenge();
};