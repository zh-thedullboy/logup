#pragma once

#include "goldilocks_quadratic_ext.h"
#include <openssl/sha.h>
#include <cstdint>
#include <vector>
#include <array>
#include <utility>

// std::array<uint8_t, 16> to_bytes(const Goldilocks2::Element& e);
// MTtree merkle_hash(const std::vector<std::vector<Goldilocks2::Element>> &data, std::array<uint8_t, SHA256_DIGEST_LENGTH> &hash);

class MerkleTree{
public:
    typedef std::array<uint8_t, SHA256_DIGEST_LENGTH> Digest;
    typedef std::vector<Digest> MTtree;
    typedef std::vector<Digest> MTPath;
    // type for a column
    typedef std::vector<Goldilocks::Element> col_t;
    typedef struct{
        MTPath path;
        col_t column;
        // this is index in the tree, not matrix!!
        size_t index;
    }MTPayload;

private:
    MTtree T;
    size_t leaf_offset;
    std::vector<col_t> cols;
public:
    MerkleTree(){};
    MerkleTree(const std::vector<std::vector<Goldilocks::Element>> &data);
    MTPayload MerkleOpen(const size_t& idx) const;
    Digest MerkleCommit() const {return T[1];}
    static bool MerkleVerify(const Digest& root, const MTPayload& payload);
};