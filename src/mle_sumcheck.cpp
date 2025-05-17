#include "mle.h"
#include "mle_sumcheck.h"
#include "goldilocks_quadratic_ext.h"
// #include <gmpxx.h>
#include <random>

sProver::sProver(const MultilinearPolynomial& g):g(g), nrnd(g.get_num_vars()), sum(Goldilocks2::zero()){
    initialize();
}

/*
1. intialize the bookkeeping table A of g in time O(l * 2 ^ l) via zeta transform
2. calculate sum
*/
void sProver::initialize() {
    uint64_t tsize = 1ull << g.get_num_vars();
    keepTable = g.get_eval_table();
    for (uint64_t mask = 0; mask < tsize; ++mask) {
        Goldilocks2::add(sum, sum, keepTable[mask]);
    }
}

std::array<Goldilocks2::Element, 2> sProver::send_message(const size_t& round, const std::vector<Goldilocks2::Element>& rands){
    std::array<Goldilocks2::Element, 2> s = {0, 0};
    
    // notation referce: Libra(https://eprint.iacr.org/2019/317.pdf) Algorithm 1
    uint64_t offset = 1ull << (nrnd - round);

    if(round > 1){
        // namely r_{i-1}
        uint64_t offset_last = (offset << 1);
        Goldilocks2::Element r = rands[round - 2];
        for(uint64_t b = 0; b < offset_last; ++b){
            Goldilocks2::Element A,B, oneminusr;
            Goldilocks2::sub(oneminusr, Goldilocks::one(), r);
            Goldilocks2::mul(A, keepTable[b], oneminusr);
            Goldilocks2::mul(B, keepTable[b + offset_last], r);
            Goldilocks2::add(keepTable[b], A, B);
        }
        keepTable.resize(offset_last);
    }
    for(uint64_t b = 0; b < offset; ++b){
        Goldilocks2::add(s[0], s[0], keepTable[b]);
        Goldilocks2::add(s[1], s[1], keepTable[b + offset]);
    }
    return s;
}

// sVerifier::sVerifier(){}

bool sVerifier::execute_sumcheck(sProver& pr, const ligeropcs_base& oracle, const size_t& sec_param){
    // if(!ligeroVerifier::check_commit(oracle, sec_param)) return false;
    Goldilocks2::Element sum = pr.get_sum();
    size_t nrnd = pr.get_rounds();
    std::vector<Goldilocks2::Element> challenges;

    // s_{i - 1}
    std::array<Goldilocks2::Element, 2> si1;
    for(int round = 1;round <= nrnd; ++round){
        // s_i
        std::array<Goldilocks2::Element, 2> si;    
        si = pr.send_message(round, challenges);
        // s(0) + s(1)
        Goldilocks2::Element ss;
        Goldilocks2::add(ss, si[0], si[1]);
        if(round == 1){
            if(!(ss == sum)) return false;
        }
        else{
            // s_{i - 1}(r) = r * s_{i - 1}(1) + (1-r) * s_{i - 1}(0)
            Goldilocks2::Element sr;
            Goldilocks2::Element r = challenges[round - 2];
            Goldilocks2::Element A, B, oneminusr;
            Goldilocks2::sub(oneminusr, Goldilocks::one(), r);
            Goldilocks2::mul(A, si1[0], oneminusr);
            Goldilocks2::mul(B, si1[1], r);
            Goldilocks2::add(sr, A, B);
            if(!(sr == ss)) return false;

            // final check
            if(round == nrnd){
                challenges.push_back(challenge());
                // should be implemented later
                // fr: f(r1, r2, ..., rl)
                Goldilocks2::Element f_r = ligeroVerifier::open(oracle, challenges, sec_param);

                // std::cout << Goldilocks2::toString(f_r) << '\n';
                // s_l(r_l)
                Goldilocks2::Element slrl;
                Goldilocks2::Element rl = challenges[round - 1];
                Goldilocks2::Element C, D, oneminusrl;
                Goldilocks2::sub(oneminusrl, Goldilocks::one(), rl);
                Goldilocks2::mul(C, si[0], oneminusrl);
                Goldilocks2::mul(D, si[1], rl);
                Goldilocks2::add(slrl, C, D);
                if(!(slrl == f_r)) return false;
                // std::cout << Goldilocks2::toString(slrl) << '\n';
            }
        }

        challenges.push_back(challenge());
        // goto next round
        si1 = si;
    }
    return true;
}

bool sVerifier::execute_sumcheck(sProver& pr, const ligeropcs_ext& oracle, const size_t& sec_param){
    // if(!ligeroVerifier::check_commit(oracle, sec_param)) return false;
    Goldilocks2::Element sum = pr.get_sum();
    size_t nrnd = pr.get_rounds();
    std::vector<Goldilocks2::Element> challenges;

    // s_{i - 1}
    std::array<Goldilocks2::Element, 2> si1;
    for(int round = 1;round <= nrnd; ++round){
        // s_i
        std::array<Goldilocks2::Element, 2> si;    
        si = pr.send_message(round, challenges);
        // s(0) + s(1)
        Goldilocks2::Element ss;
        Goldilocks2::add(ss, si[0], si[1]);
        if(round == 1){
            if(!(ss == sum)) return false;
        }
        else{
            // s_{i - 1}(r) = r * s_{i - 1}(1) + (1-r) * s_{i - 1}(0)
            Goldilocks2::Element sr;
            Goldilocks2::Element r = challenges[round - 2];
            Goldilocks2::Element A, B, oneminusr;
            Goldilocks2::sub(oneminusr, Goldilocks::one(), r);
            Goldilocks2::mul(A, si1[0], oneminusr);
            Goldilocks2::mul(B, si1[1], r);
            Goldilocks2::add(sr, A, B);
            if(!(sr == ss)) return false;

            // final check
            if(round == nrnd){
                challenges.push_back(challenge());
                // should be implemented later
                // fr: f(r1, r2, ..., rl)
                Goldilocks2::Element f_r = ligeroVerifier::open(oracle, challenges, sec_param);

                // std::cout << Goldilocks2::toString(f_r) << '\n';
                // s_l(r_l)
                Goldilocks2::Element slrl;
                Goldilocks2::Element rl = challenges[round - 1];
                Goldilocks2::Element C, D, oneminusrl;
                Goldilocks2::sub(oneminusrl, Goldilocks::one(), rl);
                Goldilocks2::mul(C, si[0], oneminusrl);
                Goldilocks2::mul(D, si[1], rl);
                Goldilocks2::add(slrl, C, D);
                if(!(slrl == f_r)) return false;
                // std::cout << Goldilocks2::toString(slrl) << '\n';
            }
        }

        challenges.push_back(challenge());
        // goto next round
        si1 = si;
    }
    return true;
}

// we use goldilocks 2-extension, so no bother specifying the field
Goldilocks2::Element sVerifier::challenge(){
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    // 设定 Goldilocks2 的模（MODULUS）为常量，假设你有这个值
    constexpr uint64_t MODULUS = Goldilocks2::p;
    std::uniform_int_distribution<uint64_t> dist(0, MODULUS - 1);

    uint64_t randn[] = {dist(gen), dist(gen)};
    return Goldilocks2::fromU64(randn);
}