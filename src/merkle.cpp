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

// hash one column
MerkleTree::Digest hash_column(const std::vector<Goldilocks2::Element>& col){
    SHA256_CTX ctx;
    SHA256_Init(&ctx);
    std::array<uint8_t, SHA256_DIGEST_LENGTH> hash;
    // 16 for 2 * 64 / 8
    for(auto e: col) SHA256_Update(&ctx, to_bytes(e).data(), 16);
    SHA256_Final(hash.data(), &ctx);
    return hash;
}

// hash two child nodes
MerkleTree::Digest hash_node(const std::array<uint8_t, SHA256_DIGEST_LENGTH>& l,const std::array<uint8_t, SHA256_DIGEST_LENGTH>& r){
    SHA256_CTX ctx;
    SHA256_Init(&ctx);
    std::array<uint8_t, SHA256_DIGEST_LENGTH> hash;
    SHA256_Update(&ctx, l.data(), SHA256_DIGEST_LENGTH);
    SHA256_Update(&ctx, r.data(), SHA256_DIGEST_LENGTH);
    SHA256_Final(hash.data(), &ctx);
    return hash;
}

// construct the merkle hash tree from a matrix
MerkleTree::MerkleTree(const std::vector<std::vector<Goldilocks2::Element>> &data){
    size_t num_cols = data[0].size();
    for(size_t i = 0; i < num_cols; ++i){
        col_t col(data.size());
        for(size_t j = 0; j < data.size(); ++j){
            col[j] = data[j][i];
        }
        cols.push_back(col);
    }

    // for clearer binary tree structure, index starts from 1 (T[0] is not used)
    T = MTtree(num_cols << 1);
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

MerkleTree::MTProof MerkleTree::MerkleOpen(const size_t& idx) const{
    // also includes the node idx itself, getting rid of the need to specify left/right
    // size_t layers = find_ceiling_log2(T.size());
    // MTPath path(layers);
    // std::array<col_t, 2> proofs_col;
    // // keep left and right
    // // 0 for left child, 1 for right child
    // proofs_col[0] = cols[idx & ~1];     // 保证 leaf_col[0] 是左节点
    // proofs_col[1] = cols[idx | 1];      // leaf_col[1] 是右节点
    size_t index = leaf_offset + idx;
    // size_t i = 0;

    MTPath path;
    while(index != 1){
        // even for left child, odd for right child
        // path[i] = T[index & ~1];
        // path[i + 1] = T[idx | 1];
        // index = index >> 1;
        // i += 2;

        size_t sibling = (index ^ 1);
        path.push_back(T[sibling]);
        index >>= 1;
    }
    return MTProof{path, cols[idx], idx + leaf_offset};
}

inline bool MerkleTree::MerkleVerify(const Digest& root, const MTProof& proof, const col_t& col){
    // Digest hashc = hash_column(col);
    // Digest hashl = hash_column(proof.leaf_col[0]);
    // Digest hashr = hash_column(proof.leaf_col[1]);
    
    // // the opened column is in the proof path
    // if(hashc != hashl && hashc != hashr) return false;

    // // verify the columns included in the proof
    // if(hashl != proof.path[0]) return false;
    // if(hashr != proof.path[1]) return false;
    
    // // verify the path
    // // idx: the index of father of the queried column in the MTTree
    // size_t i, idx = (proof.index +  (1 << (proof.path.size() >> 1))) >> 1;
    // for(i = 2; i < proof.path.size(); i += 2){
    //     if(hash_node(proof.path[i - 2], proof.path[i - 1]) != proof.path[i + (idx & 1)]) return false;
    //     idx = idx >> 2;
    // }
    // if(hash_node(proof.path[i - 2], proof.path[i - 1]) != root) return false;
    // return true;
    
    Digest hash = hash_column(col);
    size_t index = proof.index;

    for (const Digest& sibling : proof.path) {
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
    
