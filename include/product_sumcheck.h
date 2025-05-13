#pragma once

#include "goldilocks_quadratic_ext.h"
#include "mle.h"
#include "ligero.h"
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
    // typedef std::array<ligeropcs, 3> Oracle;
    static bool execute_sumcheck(pProver& pr, const std::array<ligeropcs_base, 3>& oracle, const size_t& sec_param);
    static bool execute_sumcheck(pProver& pr, const std::array<ligeropcs_ext, 3>& oracle, const size_t& sec_param);

    // customized sumcheck for \Sigma eq * frac * (gamma - p1 - lambda * p2)
    static bool execute_logup_sumcheck(
        pProver& pr,
        const MultilinearPolynomial& eqr,
        const ligeropcs_ext& frac,
        const ligeropcs_base& p1,
        const ligeropcs_base& p2,
        const Goldilocks2::Element gamma,
        const Goldilocks2::Element labmda,
        const size_t& sec_param
    );
private:
    static Goldilocks2::Element challenge();
    static inline Goldilocks2::Element mul(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e2, const  Goldilocks2::Element& e3, const  Goldilocks2::Element& e4, const  Goldilocks2::Element& e5);
    static inline Goldilocks2::Element mul(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e2, const  Goldilocks2::Element& e3);
    static inline Goldilocks2::Element add(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e2, const  Goldilocks2::Element& e3, const  Goldilocks2::Element& e4);
    static inline void interpolate_3(Goldilocks2::Element& fr, const Goldilocks2::Element& r, const Goldilocks2::Element& f1, const Goldilocks2::Element& f2, const Goldilocks2::Element& f3, const Goldilocks2::Element& f4);
};