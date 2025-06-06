#include <cstdint>
#include <string>
#include <cassert>
#include "goldilocks_quadratic_ext.h"
#include "util.h"
#include <functional>
#include <array>
#include <algorithm>
#include <iomanip>

uint64_t convert_mask_to_u64(const std::string& mask, const size_t &nvar) {
    assert(nvar >= mask.length());
    uint64_t result = 0;
    for (auto c: mask) {
        result <<= 1;
        if (c == '1') {
            result |= 1;
        }
    }
    return result;
}

int find_ceiling_log2(const uint64_t& n) {
    // at least use 1 bit
    for(int i = 1;i < 64; ++i){
        if((1ull << i) >= n){
            return i;
        }
    }
    return 64;
}

bool is_power_of_2(const uint64_t& n){
    return n != 0 && (n & (n - 1)) == 0;
}

template <typename T>
void pad_impl(std::vector<T>& table, const T& dummy){
    if(is_power_of_2(table.size())) return;
    const size_t new_size = 1ull << find_ceiling_log2(table.size());
    table.resize(new_size, dummy);
}

void pad(std::vector<Goldilocks2::Element>& table, const Goldilocks2::Element dummy) {
    pad_impl(table, dummy);
}

void pad(std::vector<Goldilocks::Element>& table, const Goldilocks::Element dummy) {
    pad_impl(table, dummy);
}

void pad(std::vector<uint64_t>& table, const uint64_t dummy) {
    pad_impl(table, dummy);
}


MultilinearPolynomial eq(const size_t& num_var, const std::vector<Goldilocks2::Element>& r){
    std::vector<Goldilocks2::Element> evaluations;
    evaluations.resize(1ull << num_var, Goldilocks2::one());
    std::vector<Goldilocks2::Element> one_minus_r(r.size());
    for (size_t i = 0; i < num_var; ++i){
        Goldilocks2::sub(one_minus_r[i], Goldilocks2::one(), r[i]);
    }
    for(uint64_t i = 0;i < num_var; ++i){
        // Goldilocks2::sub(one_minus_r[i], Goldilocks2::one(), r[i]);
        for(uint64_t j = 0;j < (1ull << i); ++j){
            // for r: low index corresponds with high bit
            Goldilocks2::mul(evaluations[j + (1ull << i)], evaluations[j], r[num_var - i - 1]);
            Goldilocks2::mul(evaluations[j], evaluations[j], one_minus_r[num_var - i - 1]);
        }
    }
    return MultilinearPolynomial(std::move(evaluations));
}

void batch_inverse(std::vector<Goldilocks2::Element>& inv, const std::vector<Goldilocks2::Element>& arr){
    assert(inv.size() >= arr.size());

    std::vector<Goldilocks2::Element> p(arr.size());
    p[0] = arr[0];
    for(size_t i = 1;i < arr.size(); ++i){
        Goldilocks2::mul(p[i], p[i - 1], arr[i]);
    }
    Goldilocks2::Element invp;
    Goldilocks2::inv(invp, p.back());
    for(size_t i = arr.size() - 1;i > 0; --i){
        Goldilocks2::mul(inv[i], invp, p[i - 1]);
        Goldilocks2::mul(invp, invp, arr[i]);
    }
    inv[0] = invp;
}

std::vector<Goldilocks::Element> random_vec_base(const size_t& n){
    srand(time(nullptr));
    std::vector<Goldilocks::Element> vec;
    vec.reserve(n);
    for(size_t i = 0;i < n; ++i){
        vec.push_back(Goldilocks::fromU64(rand()));
    }
    return vec;
}

std::vector<uint64_t> random_vec_uint(const size_t& n){
    srand(time(nullptr));
    std::vector<uint64_t> vec;
    vec.reserve(n);
    for(size_t i = 0;i < n; ++i){
        vec.push_back(rand() % RAND_MAX);
    }
    return vec;
}

std::vector<Goldilocks2::Element> random_vec_ext(const size_t& n){
    srand(time(nullptr));
    std::vector<Goldilocks2::Element> vec;
    vec.reserve(n);
    for(size_t i = 0;i < n; ++i){
        vec.push_back({Goldilocks::fromU64(rand()), Goldilocks::fromU64(rand())});
    }
    return vec;
}

std::vector<uint64_t> trange(const uint64_t& lbound, const uint64_t& ubound){
    assert(lbound <= ubound);
    std::vector<uint64_t> res;
    for(uint64_t e = lbound; e <= ubound; ++e) res.push_back(e);
    return res;
}

void alert(const std::string& mes){
    std::cout << mes << std::endl;
}

void print_table(const std::vector<Goldilocks2::Element>& table){
    for(auto e: table){
        std::cout << '(' <<  Goldilocks2::toString(e) << ')' << ' ';
    }
    std::cout << '\n';
}

void print_table(const std::vector<Goldilocks::Element>& table){
    for(auto e: table){
        std::cout << Goldilocks::toString(e) << ' ';
    }
    std::cout << '\n';
}

void print_table(const std::vector<size_t>& table){
    for(auto e: table){
        std::cout << e << ' ';
    }
    std::cout << '\n';
}

void print_hash(const std::array<uint8_t, 32>& hash){
    std::cout << "SHA256: ";
    for (int i = 0; i < 32; ++i) {
        std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
    }
    std::cout << std::endl;
}

void print_bytes(const std::array<uint8_t, 16>& hash){
    std::cout << "bytes: ";
    for (int i = 0; i < 16; ++i) {
        std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
    }
    std::cout << std::endl;
}

size_t bisearch(const std::vector<uint64_t>& arr, const uint64_t& val) {
    size_t left = 0, right = arr.size();
    while (left < right) {
        size_t mid = left + (right - left) / 2;
        if (arr[mid] < val) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }


    if (left < arr.size() && arr[left] == val) {
        return left;
    } else {
        return arr.size();
    }
}

Goldilocks2::Element Horner(const std::vector<Goldilocks2::Element> &coefs, const Goldilocks2::Element& x){
    Goldilocks2::Element res = Goldilocks2::zero();
    for(size_t i = coefs.size(); i > 0; --i){
        Goldilocks2::mul(res, res, x);
        Goldilocks2::add(res, res, coefs[i - 1]);
    }
    return res;
}

void in_place_NTT(std::vector<Goldilocks::Element>& a) {
    const size_t n = a.size();
    const size_t m = find_ceiling_log2(n);

    // Bit-reverse permutation
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    // Iterative Cooley-Tukey NTT
    for (size_t len = 2, level = 1; len <= n; len <<= 1, ++level) {
        size_t half = len >> 1;
        Goldilocks::Element w_len = ROOTS[level]; // omega^(N / len)
        for (size_t i = 0; i < n; i += len) {
            Goldilocks::Element w = Goldilocks::one();
            for (size_t j = 0; j < half; ++j) {
                auto &u = a[i + j];
                auto &v = a[i + j + half];
                Goldilocks::Element t = Goldilocks::mul(w, v);
                v = u - t;
                u = u + t;
                Goldilocks::mul(w, w, w_len);
            }
        }
    }
}

std::vector<Goldilocks::Element> NTT(const std::vector<Goldilocks::Element>& input) {
    std::vector<Goldilocks::Element> a = input;
    in_place_NTT(a);
    return a;
}

// std::vector<Goldilocks::Element> NTT(const std::vector<Goldilocks::Element>& coefs){
//     if (coefs.size() == 1) return coefs;

//     const size_t n = coefs.size();
//     const size_t m = find_ceiling_log2(n);

//     // calculate f1, f2;
//     std::vector<Goldilocks::Element> even(n >> 1), odd(n >> 1);
//     for(size_t i = 0;i < n; ++i){
//         if (i & 1) odd[i >> 1] = coefs[i];
//         else even[i >> 1] = coefs[i];
//     }
//     std::vector<Goldilocks::Element> f1 = NTT(even);
//     std::vector<Goldilocks::Element> f2 = NTT(odd);

//     std::vector<Goldilocks::Element> f(n);
//     const size_t offset = n >> 1;
//     Goldilocks::Element wn = Goldilocks::one();
//     // traverse and calculate f
//     for (size_t i = 0; i < offset; ++i){
//         f[i] = f1[i] + Goldilocks::mul(wn, f2[i]);
//         f[i + offset] = f1[i] - Goldilocks::mul(wn, f2[i]);
//         Goldilocks::mul(wn, wn, ROOTS[m]);
//     }

//     return f;
// }

size_t highest_bit_mask(const size_t& n){
    size_t mask = SIZE_MAX;
    while(mask != 0){
        if(mask & n) return mask;
        mask = mask >> 1;
    }
    return 0;
}

std::vector<Goldilocks::Element> eval_with_ntt(std::vector<Goldilocks::Element> f, const size_t& N){
    std::reverse(f.begin(), f.end());
    // N is not power of 2, pad N and trim when return
    if(!is_power_of_2(N)){
        size_t padded_N = highest_bit_mask(N) << 1;
        f.resize(padded_N, Goldilocks::zero());
        std::vector<Goldilocks::Element> raw_output = NTT(f);
        raw_output.resize(N);
        return raw_output;
    }

    // N is power of 2
    f.resize(N, Goldilocks::zero());
    return NTT(f);
}

std::vector<Goldilocks2::Element> eval_with_ntt(std::vector<Goldilocks2::Element> f, const size_t& N){
    std::reverse(f.begin(), f.end());

    size_t N_used = is_power_of_2(N) ? N : highest_bit_mask(N) << 1;
    f.resize(N_used, Goldilocks2::zero());

    std::vector<Goldilocks::Element> real(f.size()), imag(f.size());
    for(size_t i = 0; i < f.size(); ++i){
        real[i] = f[i][0];
        imag[i] = f[i][1];
    }
    
    std::vector<Goldilocks::Element> real_ntt = NTT(real);
    std::vector<Goldilocks::Element> imag_ntt = NTT(imag);

    std::vector<Goldilocks2::Element> result(N);
    for(size_t i = 0; i < N; ++i){
        result[i][0] = real_ntt[i];
        result[i][1] = imag_ntt[i];
    }
    return result;
}

std::vector<Goldilocks::Element> eval_with_ntt_base(std::vector<Goldilocks2::Element> f, const size_t& N){
    std::vector<Goldilocks::Element> base_field_copy(f.size());
    for(size_t i = 0; i < f.size(); ++i){
        base_field_copy[i] = f[i][0];
    }
    return eval_with_ntt(base_field_copy, N);
}

std::vector<Goldilocks2::Element> eval_with_ntt_ext(std::vector<Goldilocks2::Element> f, const size_t& N){
    std::vector<Goldilocks::Element> base = eval_with_ntt_base(f, N);
    std::vector<Goldilocks2::Element> ext(base.size());
    for (size_t i = 0;i < base.size(); ++i){
        ext[i][0] = base[i];
        ext[i][1] = Goldilocks::zero();
    }
    return ext;
}
