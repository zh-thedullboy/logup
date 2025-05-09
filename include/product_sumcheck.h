#pragma once

#include "goldilocks_quadratic_ext.h"
#include "mle.h"
#include <array>
#include <vector>
/*
prover for sumcheck of product three multilinear extension in O(3 * 2^l) time
*/
class pProver{
public:
    pProver(const MultilinearPolynomial& p1, const MultilinearPolynomial& p2, const MultilinearPolynomial& p3);
    void initialize();
    std::array<Goldilocks2::Element, 4> send_message(const size_t& round,const std::vector<Goldilocks2::Element>& rands);
    Goldilocks2::Element get_sum() const { return sum; }
    size_t get_rounds() const { return nrnd; }
private:
    MultilinearPolynomial p1;
    MultilinearPolynomial p2;
    MultilinearPolynomial p3;
    
    std::vector<Goldilocks2::Element> keepTablep1;
    std::vector<Goldilocks2::Element> keepTablep2;
    std::vector<Goldilocks2::Element> keepTablep3;
    inline void shrinkTable(const Goldilocks2::Element& r, const uint64_t& offset);
    static inline Goldilocks2::Element mul(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e2, const  Goldilocks2::Element& e3);
    static inline Goldilocks2::Element lincomb(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e0, const  uint64_t& r);
    Goldilocks2::Element sum;
    size_t nrnd;
};

class pVerifier{
public:
    // should be replaced with a pcs
    typedef std::array<MultilinearPolynomial, 3> Oracle;
    static bool execute_sumcheck(pProver& pr, const Oracle& oracle);
private:
    static Goldilocks2::Element challenge();
    static inline Goldilocks2::Element mul(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e2, const  Goldilocks2::Element& e3, const  Goldilocks2::Element& e4, const  Goldilocks2::Element& e5);
    static inline Goldilocks2::Element mul(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e2, const  Goldilocks2::Element& e3);
    static inline Goldilocks2::Element add(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e2, const  Goldilocks2::Element& e3, const  Goldilocks2::Element& e4);
    static inline void interpolate_3(Goldilocks2::Element& fr, const Goldilocks2::Element& r, const Goldilocks2::Element& f1, const Goldilocks2::Element& f2, const Goldilocks2::Element& f3, const Goldilocks2::Element& f4);
};