#include "product_sumcheck.h"
#include "goldilocks_quadratic_ext.h"
#include "mle.h"
#include <array>
#include <vector>
#include <random>
#include <cassert>

pProver::pProver(const MultilinearPolynomial& p1, const MultilinearPolynomial& p2, const MultilinearPolynomial& p3):p1(p1), p2(p2), p3(p3), nrnd(p1.get_num_vars()), sum(Goldilocks2::zero()){
    assert(p1.get_num_vars() == p2.get_num_vars() && p1.get_num_vars() == p3.get_num_vars());
    initialize();
}
/*
1. intialize the bookkeeping table A of g in time O(l * 2 ^ l) via zeta transform
2. calculate sum
*/
void pProver::initialize(){
    uint64_t tsize = 1ull << p1.get_num_vars();

    keepTablep1 = p1.get_eval_table();
    keepTablep2 = p2.get_eval_table();
    keepTablep3 = p3.get_eval_table();
    for (uint64_t mask = 0; mask < tsize; ++mask) {
        Goldilocks2::Element p;
        Goldilocks2::mul(p, keepTablep1[mask], keepTablep2[mask]);
        Goldilocks2::mul(p, p, keepTablep3[mask]);
        Goldilocks2::add(sum, sum, p);
    }
    // std::cout << Goldilocks2::toString(sum) << '\n';
}

inline void pProver::shrinkTable(const Goldilocks2::Element& r, const uint64_t& offset){
    for(uint64_t b = 0; b < offset; ++b){
        Goldilocks2::Element A,B, oneminusr;
        Goldilocks2::sub(oneminusr, Goldilocks::one(), r);

        Goldilocks2::mul(A, keepTablep1[b], oneminusr);
        Goldilocks2::mul(B, keepTablep1[b + offset], r);
        Goldilocks2::add(keepTablep1[b], A, B);
        
        Goldilocks2::mul(A, keepTablep2[b], oneminusr);
        Goldilocks2::mul(B, keepTablep2[b + offset], r);
        Goldilocks2::add(keepTablep2[b], A, B);

        Goldilocks2::mul(A, keepTablep3[b], oneminusr);
        Goldilocks2::mul(B, keepTablep3[b + offset], r);
        Goldilocks2::add(keepTablep3[b], A, B);
    }
    keepTablep1.resize(offset);
    keepTablep2.resize(offset);
    keepTablep3.resize(offset);
}

inline Goldilocks2::Element pProver::mul(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e2, const  Goldilocks2::Element& e3){
    Goldilocks2::Element res;
    Goldilocks2::mul(res, e1, e2);
    Goldilocks2::mul(res, res, e3);
    return res;
}


/*
linear combination
returns r * e1 + (1-r)e0
*/
inline Goldilocks2::Element pProver::lincomb(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e0, const uint64_t& r){
    Goldilocks2::Element res, A, B, oneminusr;
    Goldilocks2::sub(oneminusr, Goldilocks2::one(), r);
    Goldilocks2::mul(A, e1, r);
    Goldilocks2::mul(B, e0, oneminusr);
    Goldilocks2::add(res, A, B);
    return res;
}

std::array<Goldilocks2::Element, 4> pProver::send_message(const size_t& round, const std::vector<Goldilocks2::Element>& rands){
    std::array<Goldilocks2::Element, 4> s = {0, 0, 0, 0};
    
    // notation referce: Libra(https://eprint.iacr.org/2019/317.pdf) Algorithm 1
    uint64_t offset = 1ull << (nrnd - round);

    if(round > 1){
        // namely r_{i-1}
        uint64_t offset_last = (offset << 1);
        Goldilocks2::Element r = rands[round - 2];
        shrinkTable(r, offset_last);
    }

    for(uint64_t b = 0; b < offset; ++b){
        Goldilocks2::add(s[0], s[0], mul(keepTablep1[b], keepTablep2[b], keepTablep3[b]));
        Goldilocks2::add(s[1], s[1], mul(keepTablep1[b + offset], keepTablep2[b + offset], keepTablep3[b + offset]));
        Goldilocks2::add(s[2], s[2], mul(lincomb(keepTablep1[b + offset], keepTablep1[b], 2), lincomb(keepTablep2[b + offset], keepTablep2[b], 2), lincomb(keepTablep3[b + offset], keepTablep3[b], 2)));
        Goldilocks2::add(s[3], s[3], mul(lincomb(keepTablep1[b + offset], keepTablep1[b], 3), lincomb(keepTablep2[b + offset], keepTablep2[b], 3), lincomb(keepTablep3[b + offset], keepTablep3[b], 3)));
    }
    return s;
}

inline Goldilocks2::Element pVerifier::mul(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e2, const  Goldilocks2::Element& e3, const  Goldilocks2::Element& e4, const  Goldilocks2::Element& e5){
    Goldilocks2::Element res;
    Goldilocks2::mul(res, e1, e2);
    Goldilocks2::mul(res, res, e3);
    Goldilocks2::mul(res, res, e4);
    Goldilocks2::mul(res, res, e5);
    return res;
}
inline Goldilocks2::Element pVerifier::mul(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e2, const  Goldilocks2::Element& e3){
    Goldilocks2::Element res;
    Goldilocks2::mul(res, e1, e2);
    Goldilocks2::mul(res, res, e3);
    return res;
}
inline Goldilocks2::Element pVerifier::add(const Goldilocks2::Element& e1, const  Goldilocks2::Element& e2, const  Goldilocks2::Element& e3, const  Goldilocks2::Element& e4){
    Goldilocks2::Element res;
    Goldilocks2::add(res, e1, e2);
    Goldilocks2::add(res, res, e3);
    Goldilocks2::add(res, res, e4);
    return res;
}

/*
evaluate f(r) given f(1,2,3,4) when f is cubic
*/
inline void pVerifier::interpolate_3(Goldilocks2::Element& fr,const Goldilocks2::Element& r, const Goldilocks2::Element& f0, const Goldilocks2::Element& f1, const Goldilocks2::Element& f2, const Goldilocks2::Element& f3){
    uint64_t u6[]= {15372286724512153601ull, 0}, u2[] = {9223372034707292161ull, 0};
    Goldilocks2::Element inv6 = Goldilocks2::fromU64(u6), inv2 = Goldilocks2::fromU64(u2),minv6, minv2, x0, x1, x2, x3;
    Goldilocks2::neg(minv6, inv6);
    Goldilocks2::neg(minv2, inv2);

    Goldilocks2::sub(x0, r, 0);
    Goldilocks2::sub(x1, r, 1);
    Goldilocks2::sub(x2, r, 2);
    Goldilocks2::sub(x3, r, 3);
    fr = add(   mul(    x1, x2, x3, minv6, f0),
                mul(x0,     x2, x3, inv2,  f1),
                mul(x0, x1,     x3, minv2, f2),
                mul(x0, x1, x2,     inv6,  f3));
}

bool pVerifier::execute_sumcheck(pProver& pr, const std::array<ligeropcs_base, 3>& oracle, const size_t& sec_param){

    Goldilocks2::Element sum = pr.get_sum();
    size_t nrnd = pr.get_rounds();
    std::vector<Goldilocks2::Element> challenges;

    // s_{i - 1}
    std::array<Goldilocks2::Element, 4> si1;
    for(size_t round = 1;round <= nrnd; ++round){
        // s_i
        std::array<Goldilocks2::Element, 4> si;    
        si = pr.send_message(round, challenges);
        // s(0) + s(1)
        Goldilocks2::Element ss;
        Goldilocks2::add(ss, si[0], si[1]);
        if(round == 1){
            if(!(ss == sum)) return false;
        }
        else{
            Goldilocks2::Element sr;
            Goldilocks2::Element r = challenges[round - 2];
            interpolate_3(sr, r, si1[0], si1[1], si1[2], si1[3]);
            if(!(sr == ss)) return false;

            // final check
            if(round == nrnd){
                challenges.push_back(challenge());
                
                Goldilocks2::Element f_r = mul(ligeroVerifier::open(oracle[0], challenges, sec_param), ligeroVerifier::open(oracle[1], challenges, sec_param), ligeroVerifier::open(oracle[2], challenges, sec_param));
                Goldilocks2::Element slrl;
                Goldilocks2::Element rl = challenges[round - 1];
                interpolate_3(slrl, rl, si[0], si[1], si[2], si[3]);
                if(!(slrl == f_r)) return false;
            }
        }

        challenges.push_back(challenge());
        // goto next round
        si1 = si;
    }
    return true;
}

bool pVerifier::execute_sumcheck(pProver& pr, const std::array<ligeropcs_ext, 3>& oracle, const size_t& sec_param){

    Goldilocks2::Element sum = pr.get_sum();
    size_t nrnd = pr.get_rounds();
    std::vector<Goldilocks2::Element> challenges;

    // s_{i - 1}
    std::array<Goldilocks2::Element, 4> si1;
    for(int round = 1;round <= nrnd; ++round){
        // s_i
        std::array<Goldilocks2::Element, 4> si;    
        si = pr.send_message(round, challenges);
        // s(0) + s(1)
        Goldilocks2::Element ss;
        Goldilocks2::add(ss, si[0], si[1]);
        if(round == 1){
            if(!(ss == sum)) return false;
        }
        else{
            Goldilocks2::Element sr;
            Goldilocks2::Element r = challenges[round - 2];
            interpolate_3(sr, r, si1[0], si1[1], si1[2], si1[3]);
            if(!(sr == ss)) return false;

            // final check
            if(round == nrnd){
                challenges.push_back(challenge());
                
                Goldilocks2::Element f_r = mul(ligeroVerifier::open(oracle[0], challenges, sec_param), ligeroVerifier::open(oracle[1], challenges, sec_param), ligeroVerifier::open(oracle[2], challenges, sec_param));
                Goldilocks2::Element slrl;
                Goldilocks2::Element rl = challenges[round - 1];
                interpolate_3(slrl, rl, si[0], si[1], si[2], si[3]);
                if(!(slrl == f_r)) return false;
            }
        }

        challenges.push_back(challenge());
        // goto next round
        si1 = si;
    }
    return true;
}

bool pVerifier::execute_logup_sumcheck(
    pProver& pr,
    const MultilinearPolynomial& eqr,
    const ligeropcs_ext& frac,
    const ligeropcs_base& p1,
    const ligeropcs_base& p2,
    const Goldilocks2::Element gamma,
    const Goldilocks2::Element labmda,
    const size_t& sec_param){

    Goldilocks2::Element sum = pr.get_sum();
    size_t nrnd = pr.get_rounds();
    std::vector<Goldilocks2::Element> challenges;

    // s_{i - 1}
    std::array<Goldilocks2::Element, 4> si1;
    for(size_t round = 1;round <= nrnd; ++round){
        // s_i
        std::array<Goldilocks2::Element, 4> si;    
        si = pr.send_message(round, challenges);
        // s(0) + s(1)
        Goldilocks2::Element ss;
        Goldilocks2::add(ss, si[0], si[1]);
        if(round == 1){
            if(!(ss == sum)) return false;
        }
        else{
            Goldilocks2::Element sr;
            Goldilocks2::Element r = challenges[round - 2];
            interpolate_3(sr, r, si1[0], si1[1], si1[2], si1[3]);
            if(!(sr == ss)) return false;

            // final check, different from general product sumcheck
            if(round == nrnd){
                challenges.push_back(challenge());

                // f(r) from the oracle and the information hold by the verifier
                Goldilocks2::Element third_term;
                Goldilocks2::Element tmp;
                Goldilocks2::mul(tmp, labmda, ligeroVerifier::open(p2, challenges, sec_param));
                Goldilocks2::sub(third_term, gamma, ligeroVerifier::open(p1, challenges, sec_param));
                Goldilocks2::sub(third_term, third_term, tmp);
                Goldilocks2::Element f_r = mul(eqr.evaluate(challenges), ligeroVerifier::open(frac, challenges, sec_param), third_term);


                // f(r) from the previous rounds
                Goldilocks2::Element slrl;
                Goldilocks2::Element rl = challenges[round - 1];
                interpolate_3(slrl, rl, si[0], si[1], si[2], si[3]);
                if(!(slrl == f_r)) return false;
            }
        }

        challenges.push_back(challenge());
        // goto next round
        si1 = si;
    }
    return true;
}

// we use goldilocks 2-extension, so no bother specifying the field
Goldilocks2::Element pVerifier::challenge(){
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    
    constexpr uint64_t MODULUS = Goldilocks2::p;
    std::uniform_int_distribution<uint64_t> dist(0, MODULUS - 1);

    uint64_t randn[] = {dist(gen), dist(gen)};
    return Goldilocks2::fromU64(randn);
}