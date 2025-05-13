#include "merkle.h"
#include "../goldilocks/src/goldilocks_base_field.hpp"
#include "util.h"

#include <array>

std::array<uint8_t, 16> to_bytes(const Goldilocks2::Element& e) {
    std::array<uint8_t, 16> bytes;
    uint64_t real = Goldilocks::toU64(e[0]);
    uint64_t imag = Goldilocks::toU64(e[1]);

    for (int i = 0; i < 8; i++) {
        bytes[7 - i]  = (real >> (8 * i)) & 0xFF;
        bytes[15 - i] = (imag >> (8 * i)) & 0xFF;
    }

    // std::cout <<"bytes of ";
    // std::cout << Goldilocks2::toString(e) <<':';
    // print_bytes(bytes);
    return bytes;
}

std::array<uint8_t, 8> to_bytes(const Goldilocks::Element& e) {
    std::array<uint8_t, 8> bytes;
    uint64_t literal = Goldilocks::toU64(e);
    for (int i = 0; i < 8; i++) {
        bytes[7 - i]  = (literal >> (8 * i)) & 0xFF;
    }

    // std::cout <<"bytes of ";
    // std::cout << Goldilocks2::toString(e) <<':';
    // print_bytes(bytes);
    return bytes;
}

// hash one column
MerkleDef::Digest hash_column(const MerkleTree_ext::col_t& col){
    SHA256_CTX ctx;
    SHA256_Init(&ctx);
    std::array<uint8_t, SHA256_DIGEST_LENGTH> hash;
    // 16 for 2 * 64 / 8
    for(auto e: col) SHA256_Update(&ctx, to_bytes(e).data(), 16);
    SHA256_Final(hash.data(), &ctx);
    return hash;
}

// hash the column
MerkleDef::Digest hash_column(const MerkleTree_base::col_t& col){
    SHA256_CTX ctx;
    SHA256_Init(&ctx);
    std::array<uint8_t, SHA256_DIGEST_LENGTH> hash;
    // 8 for 64 / 8
    for(auto e: col) SHA256_Update(&ctx, to_bytes(e).data(), 8);
    SHA256_Final(hash.data(), &ctx);
    return hash;
}

// hash two child nodes
MerkleDef::Digest hash_node(const std::array<uint8_t, SHA256_DIGEST_LENGTH>& l,const std::array<uint8_t, SHA256_DIGEST_LENGTH>& r){
    SHA256_CTX ctx;
    SHA256_Init(&ctx);
    std::array<uint8_t, SHA256_DIGEST_LENGTH> hash;
    SHA256_Update(&ctx, l.data(), SHA256_DIGEST_LENGTH);
    SHA256_Update(&ctx, r.data(), SHA256_DIGEST_LENGTH);
    SHA256_Final(hash.data(), &ctx);
    return hash;
}

// construct the merkle hash tree from a matrix
MerkleTree_base::MerkleTree_base(const std::vector<col_t> &data){
    size_t num_cols = data[0].size();
    for(size_t i = 0; i < num_cols; ++i){
        col_t col(data.size());
        for(size_t j = 0; j < data.size(); ++j){
            col[j] = data[j][i];
        }
        cols.push_back(col);
    }

    // for clearer binary tree structure, index starts from 1 (T[0] is not used)
    T = MerkleDef::MTtree(num_cols << 1);
    // MTtree mt_t(num_cols << 1);
    size_t loopN = find_ceiling_log2(num_cols);
    leaf_offset = 1ul << loopN;
    for(size_t j = 0; j < leaf_offset; ++j){
        T[leaf_offset + j] = hash_column(cols[j]);
    }
    for(size_t i = loopN; i > 0; --i){
        size_t offset = 1ul << (i - 1);
        for(size_t j = 0; j < offset; ++j){
            size_t idx = offset + j;
            T[idx] = hash_node(T[2 * idx], T[2 * idx + 1]);
        }
    }
}

MerkleTree_base::MTPayload MerkleTree_base::MerkleOpen(const size_t& idx) const{
    size_t index = leaf_offset + idx;
    MerkleDef::MTPath path;
    while(index != 1){
        size_t sibling = (index ^ 1);
        path.push_back(T[sibling]);
        index >>= 1;
    }
    return MTPayload{path, cols[idx], idx + leaf_offset};
}

bool MerkleTree_base::MerkleVerify(const MerkleDef::Digest& root, const MTPayload& payload){
    col_t col = payload.column;
    MerkleDef::Digest hash = hash_column(col);
    size_t index = payload.index;

    for (const MerkleDef::Digest& sibling : payload.path) {
        if ((index & 1) == 0) {
            // this is left child
            hash = hash_node(hash, sibling);
        } else {
            // this is right child
            hash = hash_node(sibling, hash);
        }
        index >>= 1;
    }
    return hash == root;
}


// construct the merkle hash tree from a matrix
MerkleTree_ext::MerkleTree_ext(const std::vector<col_t> &data){
    size_t num_cols = data[0].size();
    for(size_t i = 0; i < num_cols; ++i){
        col_t col(data.size());
        for(size_t j = 0; j < data.size(); ++j){
            col[j] = data[j][i];
        }
        cols.push_back(col);
    }

    // for clearer binary tree structure, index starts from 1 (T[0] is not used)
    T = MerkleDef::MTtree(num_cols << 1);
    // MTtree mt_t(num_cols << 1);
    size_t loopN = find_ceiling_log2(num_cols);
    leaf_offset = 1ul << loopN;
    for(size_t j = 0; j < leaf_offset; ++j){
        T[leaf_offset + j] = hash_column(cols[j]);
    }
    for(size_t i = loopN; i > 0; --i){
        size_t offset = 1ul << (i - 1);
        for(size_t j = 0; j < offset; ++j){
            size_t idx = offset + j;
            T[idx] = hash_node(T[2 * idx], T[2 * idx + 1]);
        }
    }
}

MerkleTree_ext::MTPayload MerkleTree_ext::MerkleOpen(const size_t& idx) const{
    size_t index = leaf_offset + idx;
    MerkleDef::MTPath path;
    while(index != 1){
        size_t sibling = (index ^ 1);
        path.push_back(T[sibling]);
        index >>= 1;
    }
    return MTPayload{path, cols[idx], idx + leaf_offset};
}

bool MerkleTree_ext::MerkleVerify(const MerkleDef::Digest& root, const MTPayload& payload){
    col_t col = payload.column;
    MerkleDef::Digest hash = hash_column(col);
    size_t index = payload.index;

    for (const MerkleDef::Digest& sibling : payload.path) {
        if ((index & 1) == 0) {
            // this is left child
            hash = hash_node(hash, sibling);
        } else {
            // this is right child
            hash = hash_node(sibling, hash);
        }
        index >>= 1;
    }
    return hash == root;
}
