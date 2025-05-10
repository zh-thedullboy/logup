#include "ligero.h"
#include "mle.h"
#include "goldilocks_quadratic_ext.h"
#include "merkle.h"
#include "util.h"
#include <cmath>
#include <openssl/sha.h>
#include <cassert>

std::vector<Goldilocks2::Element> rsencode(const std::vector<Goldilocks2::Element> &data, const double& rho){
    size_t n = std::ceil(data.size() / rho);
    std::vector<Goldilocks2::Element> code;
    code.reserve(n);
    for(size_t i = 0; i < n; ++i){
        code.push_back(Horner(data, Goldilocks2::fromU64(i)));
    }
    return code;
}

ligeroProver::ligeroProver(const MultilinearPolynomial& w, const double& rho):mle(w), rho(rho){
    // stevals = w.get_eval_table();
    std::vector<Goldilocks2::Element> evals = w.get_eval_table();
    size_t l = w.get_num_vars();

    // 2^l = a * b
    a = 1ull << (l >> 1);       //floor(l/2)
    b = a << (l & 1);           //ceil(l/2)
    codelen = std::ceil(b / rho);
    std::vector<std::vector<Goldilocks2::Element>> beforerscode;
    for(size_t i = 0; i < a; ++i){
        std::vector<Goldilocks2::Element> dataline(b);
        for(size_t j = 0; j < b; ++j){
            // dataline.push_back(evals[i * b + j]);
            dataline[j] = evals[i * b + j];
        }
        beforerscode.push_back(dataline);
        codewords.push_back(rsencode(dataline, rho));
    }
    mt_t = MerkleTree(codewords);
}

std::vector<Goldilocks2::Element> ligeroProver::lincomb(const std::vector<Goldilocks2::Element>& r) const{
    assert(r.size() == a);
    // std::cout << r.size() << '\n' << a << '\n';
    std::vector<Goldilocks2::Element> v(b, Goldilocks2::zero());
    std::vector<Goldilocks2::Element> evals = mle.get_eval_table();
    for(size_t j = 0; j < a; ++j){
        Goldilocks2::Element tmp;
        size_t offset = j * b;
        for(size_t i = 0; i < b; ++i){
            Goldilocks2::mul(tmp, r[j], evals[offset + i]);
            Goldilocks2::add(v[i], v[i], tmp);
        }
    }
    return v;
}

std::vector<MerkleTree::MTPayload> ligeroProver::open_cols(const std::vector<size_t>& indexes) const{
    std::vector<MerkleTree::MTPayload> payloads;
    for(size_t e: indexes){
        payloads.push_back(mt_t.MerkleOpen(e));
    }
    return payloads;
}

ligeropcs ligeroProver::commit() const{
    return {mt_t.MerkleCommit(), std::make_shared<ligeroProver>(*this), a, b};
}

std::mt19937_64 ligeroVerifier::gen(std::random_device{}());
std::uniform_int_distribution<uint64_t> ligeroVerifier::dist(0, Goldilocks2::p - 1);

Goldilocks2::Element ligeroVerifier::randnum(){
    return {Goldilocks::fromU64(dist(gen)), Goldilocks::fromU64(dist(gen))};
}

std::vector<Goldilocks2::Element> ligeroVerifier::randvec(const uint64_t& n){
    std::vector<Goldilocks2::Element> rands;
    rands.reserve(n);
    for (uint64_t i = 0; i < n; ++i) {
        rands.push_back(randnum());
    }
    return rands;
}

std::vector<size_t> ligeroVerifier::randindexes(const uint64_t& n, const size_t& bound){
    static std::mt19937 _gen(std::random_device{}());
    std::uniform_int_distribution<size_t> _dist(0, bound - 1);
    std::vector<size_t> rands;
    rands.reserve(n);
    for (uint64_t i = 0; i < n; ++i) {
        rands.push_back(_dist(_gen));
    }
    return rands;
}

bool ligeroVerifier::check_lincomb(const ligeropcs& pcs, const std::vector<Goldilocks2::Element>& r , const std::vector<Goldilocks2::Element>& comb){
    const ligeroProver &prover = *pcs.prover;
    // 要 t 个 column   !!!t要改成 \theta(\lambda)
    size_t t = 5;
    std::vector<size_t> indexes = randindexes(t, std::ceil(pcs.num_cols / prover.rho));
    // std::vector<size_t> indexes = {7, 7, 7, 7, 7};

    std::vector<MerkleTree::MTPayload> openings = prover.open_cols(indexes);
    for(auto opening: openings){
        // check if this opening is right
        if(!MerkleTree::MerkleVerify(pcs.mthash, opening)) std::cout <<"got ya\n";
    }

    for(size_t k = 0;k < indexes.size(); ++k){
        size_t idx = indexes[k];
        // check if this entry is correctly computed
        Goldilocks2::Element entry = Goldilocks2::zero();
        Goldilocks2::Element tmp;
        for(size_t i = 0;i < pcs.num_rows; ++i){
            Goldilocks2::mul(tmp, r[i], openings[k].column[i]);
            Goldilocks2::add(entry, entry, tmp);
        }
        if(entry != comb[idx]) return false;
    }

    return true;
}

bool ligeroVerifier::check_commit(const ligeropcs& pcs){
    const ligeroProver& prover = *pcs.prover; 
    std::vector<Goldilocks2::Element> r = randvec(pcs.num_rows);
    std::vector<Goldilocks2::Element> v = prover.lincomb(r);
    std::vector<Goldilocks2::Element> w = rsencode(v, prover.rho);

    return check_lincomb(pcs, r, w);
}

// maybe can be moved to util?
// calculate b^T dot u dot a, notation from pazk(https://people.cs.georgetown.edu/jthaler/ProofsArgsAndZK.pdf)
Goldilocks2::Element tensor_product(const std::vector<Goldilocks2::Element> &b, const std::vector<std::vector<Goldilocks2::Element>>& u, const std::vector<Goldilocks2::Element> &a){
    // assume (1, h) * (h, w) * (w, 1)
    size_t h = b.size();
    size_t w = a.size();
    assert(h == u.size() && u[0].size() == w);
    Goldilocks2::Element res = Goldilocks2::zero();
    for(size_t i = 0; i < w; ++i){
        Goldilocks2::Element term = Goldilocks2::zero();
        Goldilocks2::Element tmp;
        for(size_t j = 0; j < h; ++j){
            Goldilocks2::mul(tmp, b[j], u[j][i]);
            Goldilocks2::add(term, term, tmp);
        }
        Goldilocks2::add(res, res, term);
    }
    return res;
}

Goldilocks2::Element dot_product(const std::vector<Goldilocks2::Element> &b, const std::vector<Goldilocks2::Element> &a){
    size_t n = b.size();
    assert(n == a.size());
    Goldilocks2::Element res = Goldilocks2::zero();
    for(size_t i = 0; i < n; ++i){
        Goldilocks2::Element tmp;
        Goldilocks2::mul(tmp, b[i], a[i]);
        Goldilocks2::add(res, res, tmp);
    }
    return res;
}

Goldilocks2::Element ligeroVerifier::open(const ligeropcs& pcs, const std::vector<Goldilocks2::Element> &z){
    // 自己算L和R
    std::array<std::vector<Goldilocks2::Element>, 2> lr = calculate_lr(z.size(), z);
    std::vector<Goldilocks2::Element> L = lr[0];
    std::vector<Goldilocks2::Element> R = lr[1];

    const ligeroProver &prover = *pcs.prover;
    // v_prime for v', idealy we have E(v') = w', w' = R dot uhat
    std::vector<Goldilocks2::Element> v_prime = prover.lincomb(R);

    std::vector<Goldilocks2::Element> w_prime = rsencode(v_prime, prover.rho);
    assert(check_lincomb(pcs, R, w_prime));

    return dot_product(v_prime, L);
}

std::array<std::vector<Goldilocks2::Element>, 2> ligeroVerifier::calculate_lr(const size_t& num_var, const std::vector<Goldilocks2::Element> &z){
    //different from a,b in prover, a, b here are the log of each
    size_t a = num_var >> 1, b = a + (num_var & 1);

    // high bit of z
    std::vector<Goldilocks2::Element> zh(z.begin(), z.begin() + a);
    // low bits of z
    std::vector<Goldilocks2::Element> zl(z.begin() + a, z.end());

    std::vector<Goldilocks2::Element> L = eq(b, zl).get_eval_table();
    std::vector<Goldilocks2::Element> R = eq(a, zh).get_eval_table();

    return {L, R};
}

