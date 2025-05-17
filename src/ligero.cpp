#include "ligero.h"
// #include "mle.h"
#include "goldilocks_quadratic_ext.h"
#include "merkle.h"
#include "util.h"
#include "timer.h"
#include <cmath>
#include <openssl/sha.h>
#include <cassert>

// reed solomon encode data on base field
std::vector<Goldilocks::Element> rsencode(const std::vector<Goldilocks::Element> &data, const uint64_t& rho_inv){
    // return eval_with_ntt_ext(data, data.size() * rho_inv);
    return eval_with_ntt(data, data.size() * rho_inv);
}

// reed solomon encode data on quadratic extension field
// this is simply encode real part and imaginary part respectively and merge the results
std::vector<Goldilocks2::Element> rsencode(const std::vector<Goldilocks2::Element> &data, const uint64_t& rho_inv){
    return eval_with_ntt(data, data.size() * rho_inv);
}

ligeroProver_base::ligeroProver_base(const MultilinearPolynomial& w, const uint64_t& rho_inv):rho_inv(rho_inv){
    size_t l = w.get_num_vars();

    // 2^l = a * b
    a = 1ull << (l >> 1);       //floor(l/2)
    b = a << (l & 1);           //ceil(l/2)
    M.resize(1ull << l, Goldilocks::zero());
    codelen = b * rho_inv;
    for(size_t i = 0; i < w.get_eval_table().size(); ++i){
        M[i] = w.eval_hypercube(i)[0];
    }

    for(size_t i = 0; i < a; ++i){
        std::vector<Goldilocks::Element> dataline(b);
        for(size_t j = 0; j < b; ++j){
            dataline[j] = M[i * b + j];
        }
        codewords.push_back(rsencode(dataline, rho_inv));
    }
    mt_t = MerkleTree_base(codewords);
}


ligeroProver_base::ligeroProver_base(const std::vector<Goldilocks::Element>& w, const uint64_t& rho_inv):M(w), rho_inv(rho_inv){
    // stevals = w.get_eval_table();
    size_t l = find_ceiling_log2(w.size());

    // 2^l = a * b
    a = 1ull << (l >> 1);       //floor(l/2)
    b = a << (l & 1);           //ceil(l/2)
    M.resize(1ull << l, Goldilocks::zero());
    codelen = b * rho_inv;
    for(size_t i = 0; i < w.size(); ++i){
        M[i] = w[i];
    }
    for(size_t i = 0; i < a; ++i){
        std::vector<Goldilocks::Element> dataline(b);
        for(size_t j = 0; j < b; ++j){
            dataline[j] = M[i * b + j];
        }
        codewords.push_back(rsencode(dataline, rho_inv));
    }
    mt_t = MerkleTree_base(codewords);
}

ligeroProver_base::ligeroProver_base(const std::vector<uint64_t>& w, const uint64_t& rho_inv):rho_inv(rho_inv){
    // stevals = w.get_eval_table();
    size_t l = find_ceiling_log2(w.size());
    // std::cout << l << '\n';

    // 2^l = a * b
    a = 1ull << (l >> 1);       //floor(l/2)
    b = a << (l & 1);           //ceil(l/2)
    set_timer("make matrix");
    M.resize(1ull << l, Goldilocks::zero());
    codelen = b * rho_inv;
    for(size_t i = 0; i < w.size(); ++i){
        M[i] = Goldilocks::fromU64(w[i]);
    }
    set_timer("ntt");
    pause_timer("ntt");
    alert("log-size of vector to be ntt-ed: " + std::to_string(b));
    for(size_t i = 0; i < a; ++i){
        std::vector<Goldilocks::Element> dataline(b);
        for(size_t j = 0; j < b; ++j){
            dataline[j] = M[i * b + j];
        }
        resume_timer("ntt");
        codewords.push_back(rsencode(dataline, rho_inv));
        pause_timer("ntt");
    }
    end_timer("ntt");
    end_timer("make matrix");
    set_timer("build merkle tree");
    mt_t = MerkleTree_base(codewords);
    end_timer("build merkle tree");
}

std::vector<Goldilocks2::Element> ligeroProver_base::lincomb(const std::vector<Goldilocks2::Element>& r) const{
    assert(r.size() == a);
    // std::cout << r.size() << '\n' << a << '\n';
    std::vector<Goldilocks2::Element> v(b, Goldilocks2::zero());
    for(size_t j = 0; j < a; ++j){
        Goldilocks2::Element tmp;
        size_t offset = j * b;
        for(size_t i = 0; i < b; ++i){
            Goldilocks2::mul(tmp, r[j], M[offset + i]);
            Goldilocks2::add(v[i], v[i], tmp);
        }
    }
    return v;
}

std::vector<MerkleTree_base::MTPayload> ligeroProver_base::open_cols(const std::vector<size_t>& indexes) const{
    std::vector<MerkleTree_base::MTPayload> payloads;
    for(size_t e: indexes){
        payloads.push_back(mt_t.MerkleOpen(e));
    }
    return payloads;
}

ligeropcs_base ligeroProver_base::commit() const{
    return {mt_t.MerkleCommit(), std::make_shared<ligeroProver_base>(*this), a, b};
}



ligeroProver_ext::ligeroProver_ext(const MultilinearPolynomial& w, const uint64_t& rho_inv):rho_inv(rho_inv){
    // evals = w.get_eval_table();
    // std::vector<Goldilocks2::Element> evals = w.get_eval_table();
    size_t l = w.get_num_vars();
    M.resize(1ull << l, Goldilocks2::zero());
    // 2^l = a * b
    a = 1ull << (l >> 1);       //floor(l/2)
    b = a << (l & 1);           //ceil(l/2)
    codelen = b * rho_inv;
    for(size_t i = 0; i < w.get_eval_table().size(); ++i){
        M[i] = w.eval_hypercube(i);
    }
    for(size_t i = 0; i < a; ++i){
        std::vector<Goldilocks2::Element> dataline(b);
        for(size_t j = 0; j < b; ++j){
            dataline[j] = M[i * b + j];
        }
        codewords.push_back(rsencode(dataline, rho_inv));
    }
    mt_t = MerkleTree_ext(codewords);
}

ligeroProver_ext::ligeroProver_ext(const std::vector<Goldilocks2::Element>& w, const uint64_t& rho_inv):M(w), rho_inv(rho_inv){
    size_t l = find_ceiling_log2(w.size());
    M.resize(1ull << l, Goldilocks2::zero());
    // 2^l = a * b
    a = 1ull << (l >> 1);       //floor(l/2)
    b = a << (l & 1);           //ceil(l/2)
    codelen = b * rho_inv;
    for(size_t i = 0; i < w.size(); ++i){
        M[i] = w[i];
    }
    for(size_t i = 0; i < a; ++i){
        std::vector<Goldilocks2::Element> dataline(b);
        for(size_t j = 0; j < b; ++j){
            dataline[j] = M[i * b + j];
        }
        codewords.push_back(rsencode(dataline, rho_inv));
    }
    mt_t = MerkleTree_ext(codewords);
}

std::vector<Goldilocks2::Element> ligeroProver_ext::lincomb(const std::vector<Goldilocks2::Element>& r) const{
    assert(r.size() == a);
    // std::cout << r.size() << '\n' << a << '\n';
    std::vector<Goldilocks2::Element> v(b, Goldilocks2::zero());
    for(size_t j = 0; j < a; ++j){
        Goldilocks2::Element tmp;
        size_t offset = j * b;
        for(size_t i = 0; i < b; ++i){
            Goldilocks2::mul(tmp, r[j], M[offset + i]);
            Goldilocks2::add(v[i], v[i], tmp);
        }
    }
    return v;
}

std::vector<MerkleTree_ext::MTPayload> ligeroProver_ext::open_cols(const std::vector<size_t>& indexes) const{
    std::vector<MerkleTree_ext::MTPayload> payloads;
    for(size_t e: indexes){
        payloads.push_back(mt_t.MerkleOpen(e));
    }
    return payloads;
}

ligeropcs_ext ligeroProver_ext::commit() const{
    return {mt_t.MerkleCommit(), std::make_shared<ligeroProver_ext>(*this), a, b};
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


bool ligeroVerifier::check_lincomb(const ligeropcs_base& pcs, const std::vector<Goldilocks2::Element>& r , const std::vector<Goldilocks2::Element>& comb, const size_t& t){
    const auto &prover = *pcs.prover;
    std::vector<size_t> indexes = randindexes(t, std::ceil(pcs.num_cols * prover.rho_inv));
    // std::vector<size_t> indexes = {7, 7, 7, 7, 7};

    auto openings = prover.open_cols(indexes);
    for(auto opening: openings){
        // check if this opening is right
        if(!MerkleTree_base::MerkleVerify(pcs.mthash, opening)) return false;
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

bool ligeroVerifier::check_lincomb(const ligeropcs_ext& pcs, const std::vector<Goldilocks2::Element>& r , const std::vector<Goldilocks2::Element>& comb, const size_t& t){
    const auto &prover = *pcs.prover;
    std::vector<size_t> indexes = randindexes(t, pcs.num_cols * prover.rho_inv);
    // std::vector<size_t> indexes = {7, 7, 7, 7, 7};

    auto openings = prover.open_cols(indexes);
    for(auto opening: openings){
        // check if this opening is right
        if(!MerkleTree_ext::MerkleVerify(pcs.mthash, opening)) return false;
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

bool ligeroVerifier::check_commit(const ligeropcs_base& pcs, const size_t& sec_param){
    const auto& prover = *pcs.prover; 
    std::vector<Goldilocks2::Element> r = randvec(pcs.num_rows);
    std::vector<Goldilocks2::Element> v = prover.lincomb(r);
    std::vector<Goldilocks2::Element> w = rsencode(v, prover.rho_inv);
    size_t t = calculate_t(sec_param, prover.rho_inv, prover.codelen, FIELD_BITS);
    return check_lincomb(pcs, r, w, t);
}

bool ligeroVerifier::check_commit(const ligeropcs_ext& pcs, const size_t& sec_param){
    const auto& prover = *pcs.prover; 
    std::vector<Goldilocks2::Element> r = randvec(pcs.num_rows);
    std::vector<Goldilocks2::Element> v = prover.lincomb(r);
    std::vector<Goldilocks2::Element> w = rsencode(v, prover.rho_inv);
    size_t t = calculate_t(sec_param, prover.rho_inv, prover.codelen, FIELD_BITS);
    return check_lincomb(pcs, r, w, t);
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

Goldilocks2::Element ligeroVerifier::open(const ligeropcs_base& pcs, const std::vector<Goldilocks2::Element> &z,  const size_t& sec_param){
    std::array<std::vector<Goldilocks2::Element>, 2> lr = calculate_lr(z.size(), z);
    
    std::vector<Goldilocks2::Element> L = lr[0];
    std::vector<Goldilocks2::Element> R = lr[1];

    const auto &prover = *pcs.prover;
    // v_prime for v', idealy we have E(v') = w', w' = R dot uhat
    std::vector<Goldilocks2::Element> v_prime = prover.lincomb(R);

    std::vector<Goldilocks2::Element> w_prime = rsencode(v_prime, prover.rho_inv);
    size_t t = calculate_t(sec_param, prover.rho_inv, prover.codelen, FIELD_BITS);
    assert(check_lincomb(pcs, R, w_prime, t));

    return dot_product(v_prime, L);
}

Goldilocks2::Element ligeroVerifier::open(const ligeropcs_ext& pcs, const std::vector<Goldilocks2::Element> &z,  const size_t& sec_param){
    std::array<std::vector<Goldilocks2::Element>, 2> lr = calculate_lr(z.size(), z);
    
    std::vector<Goldilocks2::Element> L = lr[0];
    std::vector<Goldilocks2::Element> R = lr[1];

    const auto &prover = *pcs.prover;
    // v_prime for v', idealy we have E(v') = w', w' = R dot uhat
    std::vector<Goldilocks2::Element> v_prime = prover.lincomb(R);

    std::vector<Goldilocks2::Element> w_prime = rsencode(v_prime, prover.rho_inv);
    size_t t = calculate_t(sec_param, prover.rho_inv, prover.codelen, FIELD_BITS);
    assert(check_lincomb(pcs, R, w_prime, t));

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

size_t ligeroVerifier::calculate_t(
    const size_t& sec_param,
    const uint64_t& rho_inv,
    const size_t& codeword_len,
    const size_t& field_bits
) {
    // residual = n / 2^field_bits
    double residual = static_cast<double>(codeword_len) / std::pow(2.0, field_bits);

    double pr = std::pow(2.0, -static_cast<int>(sec_param));
    // unless target sec level can't be achieved
    assert(pr > residual);

    double numerator = std::log2(pr - residual) - 1;
    double denominator = std::log2((1.0 + 1.0 / static_cast<double>(rho_inv)) * 0.5);
    size_t t = static_cast<size_t>(std::ceil(numerator / denominator));
    return std::min(t, codeword_len);
}
