#pragma once

#include "mle.h"
#include <cstdint>
#include <vector>
#include "merkle.h"
#include <memory>
#include <random>

#define FIELD_BITS 2 * 64
class ligeroProver_base;
class ligeroProver_ext;

typedef struct{
    MerkleDef::Digest mthash;
    // const ligeroProver& prover;
    std::shared_ptr<ligeroProver_base> prover;
    size_t num_rows;
    size_t num_cols;
}ligeropcs_base;

typedef struct{
    MerkleDef::Digest mthash;
    // const ligeroProver& prover;
    std::shared_ptr<ligeroProver_ext> prover;
    size_t num_rows;
    size_t num_cols;
}ligeropcs_ext;

class ligeroProver_base{
public:
    // code rate
    uint64_t rho_inv;
    size_t codelen;

public:
    ligeroProver_base(const MultilinearPolynomial& w, const uint64_t& rho_inv);
    ligeroProver_base(const std::vector<uint64_t>& w, const uint64_t& rho_inv);
    ligeroProver_base(const std::vector<Goldilocks::Element>& w, const uint64_t& rho_inv);
    ligeropcs_base commit() const;
    std::vector<Goldilocks2::Element> lincomb(const std::vector<Goldilocks2::Element>& r) const;
    std::vector<MerkleTree_base::MTPayload> open_cols(const std::vector<size_t>& indexes) const;
    
private:
    // mle is a multilinear polynomial whose evaluations over the hypercube are all base field elements 
    // MultilinearPolynomial mle;

    // num of rows, columns;
    size_t a, b;
    // original vector
    std::vector<Goldilocks::Element> M;
    // encoded matrix
    std::vector<std::vector<Goldilocks::Element>> codewords;
    // merkle hash tree of codewords
    MerkleTree_base mt_t;
};

class ligeroProver_ext{
public:
    // code rate
    uint64_t rho_inv;
    size_t codelen;

public:
    ligeroProver_ext(const MultilinearPolynomial& w, const uint64_t& rho_inv);
    ligeroProver_ext(const std::vector<Goldilocks2::Element>& w, const uint64_t& rho_inv);
    ligeropcs_ext commit() const;
    std::vector<Goldilocks2::Element> lincomb(const std::vector<Goldilocks2::Element>& r) const;
    std::vector<MerkleTree_ext::MTPayload> open_cols(const std::vector<size_t>& indexes) const;
    
private:
    // MultilinearPolynomial mle;

    // num of rows, columns;
    size_t a, b;
    // original vector
    std::vector<Goldilocks2::Element> M;
    // encoded matrix
    std::vector<std::vector<Goldilocks2::Element>> codewords;
    // merkle hash tree of codewords
    MerkleTree_ext mt_t;
};



class ligeroVerifier{
public:
    // check if some commit is valid ligero commit
    
    static bool check_commit(const ligeropcs_base& pcs, const size_t& sec_param);
    static bool check_commit(const ligeropcs_ext& pcs, const size_t& sec_param);
    // static bool check_commit(const ligeropcs& pcs, const size_t& sec_param);

    // open f(z) where f is a polynomial hold by prover
    // verification is included in this process
    static Goldilocks2::Element open(const ligeropcs_base& pcs, const std::vector<Goldilocks2::Element> &z,  const size_t& sec_param);
    static Goldilocks2::Element open(const ligeropcs_ext& pcs, const std::vector<Goldilocks2::Element> &z,  const size_t& sec_param);
private:
    static std::mt19937_64 gen;
    static std::uniform_int_distribution<uint64_t> dist;
    static Goldilocks2::Element randnum();
    static std::vector<Goldilocks2::Element> randvec(const uint64_t& n);
    static std::vector<size_t> randindexes(const uint64_t& n, const size_t& bound);
    static std::array<std::vector<Goldilocks2::Element>, 2> calculate_lr(const size_t& num_var, const std::vector<Goldilocks2::Element> &z);
    static bool check_lincomb(const ligeropcs_base& pcs, const std::vector<Goldilocks2::Element>& r, const std::vector<Goldilocks2::Element>& comb, const size_t& t);
    static bool check_lincomb(const ligeropcs_ext& pcs, const std::vector<Goldilocks2::Element>& r, const std::vector<Goldilocks2::Element>& comb, const size_t& t);
    static size_t calculate_t(const size_t& sec_param, const uint64_t& rho_inv, const size_t& codeword_len, const size_t& field_bits);
};



std::vector<Goldilocks2::Element> rsencode(const std::vector<Goldilocks2::Element> &data, const uint64_t& rho_inv);
std::vector<Goldilocks::Element> rsencode(const std::vector<Goldilocks::Element> &data, const uint64_t& rho_inv);