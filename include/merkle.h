#pragma once

#include "goldilocks_quadratic_ext.h"
#include <openssl/sha.h>
#include <cstdint>
#include <vector>
#include <array>
#include <utility>

// std::array<uint8_t, 16> to_bytes(const Goldilocks2::Element& e);
// MTtree merkle_hash(const std::vector<std::vector<Goldilocks2::Element>> &data, std::array<uint8_t, SHA256_DIGEST_LENGTH> &hash);

namespace MerkleDef{
    typedef std::array<uint8_t, SHA256_DIGEST_LENGTH> Digest;
    typedef std::vector<Digest> MTtree;
    typedef std::vector<Digest> MTPath;
}

class MerkleTree_base{
public:

    // type for a column
    typedef std::vector<Goldilocks::Element> col_t;
    typedef struct{
        MerkleDef::MTPath path;
        col_t column;
        // this is index in the tree, not matrix!!
        size_t index;
    }MTPayload;

private:
    MerkleDef::MTtree T;
    size_t leaf_offset;
    std::vector<col_t> cols;
public:
    MerkleTree_base(){};
    MerkleTree_base(const std::vector<std::vector<Goldilocks::Element>> &data);
    MTPayload MerkleOpen(const size_t& idx) const;
    MerkleDef::Digest MerkleCommit() const {return T[1];}
    static bool MerkleVerify(const MerkleDef::Digest& root, const MTPayload& payload);
};

class MerkleTree_ext{
public:
    // type for a column
    typedef std::vector<Goldilocks2::Element> col_t;
    typedef struct{
        MerkleDef::MTPath path;
        col_t column;
        // this is index in the tree, not matrix!!
        size_t index;
    }MTPayload;

private:
    MerkleDef::MTtree T;
    size_t leaf_offset;
    std::vector<col_t> cols;
public:
    MerkleTree_ext(){};
    MerkleTree_ext(const std::vector<std::vector<Goldilocks2::Element>> &data);
    MTPayload MerkleOpen(const size_t& idx) const;
    MerkleDef::Digest MerkleCommit() const {return T[1];}
    static bool MerkleVerify(const MerkleDef::Digest& root, const MTPayload& payload);
};