#pragma once

#include "mle.h"
#include <cstdint>
#include <vector>
#include "merkle.h"
#include <memory>
#include <random>

class ligeroProver;

typedef struct{
    MerkleTree::Digest mthash;
    // const ligeroProver& prover;
    std::shared_ptr<ligeroProver> prover;
    size_t num_rows;
    size_t num_cols;
}ligeropcs;

class ligeroProver{
public:
    // code rate
    double rho;

public:
    ligeroProver(const MultilinearPolynomial& w, const double& rho);
    ligeropcs commit() const;
    std::vector<Goldilocks2::Element> lincomb(const std::vector<Goldilocks2::Element>& r) const;
    std::vector<MerkleTree::MTPayload> open_cols(const std::vector<size_t>& indexes) const; 
    
private:
    MultilinearPolynomial mle;

    // num of rows, columns and length of a code;
    size_t a, b, codelen;
    // encoded matrix
    std::vector<std::vector<Goldilocks2::Element>> codewords;
    // merkle hash tree of codewords
    MerkleTree mt_t;
};

class ligeroVerifier{
public:
    // check if some commit is valid ligero commit
    static bool check_commit(const ligeropcs& pcs);

    // open f(z) where f is a polynomial hold by prover
    // verification is included in this process
    static Goldilocks2::Element open(const ligeropcs& pcs, const std::vector<Goldilocks2::Element> &z);
private:
    static std::mt19937_64 gen;
    static std::uniform_int_distribution<uint64_t> dist;
    static Goldilocks2::Element randnum();
    static std::vector<Goldilocks2::Element> randvec(const uint64_t& n);
    static std::vector<size_t> randindexes(const uint64_t& n, const size_t& bound);
    static std::array<std::vector<Goldilocks2::Element>, 2> calculate_lr(const size_t& num_var, const std::vector<Goldilocks2::Element> &z);
    static bool check_lincomb(const ligeropcs& pcs, const std::vector<Goldilocks2::Element>& r , const std::vector<Goldilocks2::Element>& comb);
};

std::vector<Goldilocks2::Element> rsencode(const std::vector<Goldilocks2::Element> &data, const double& rho);