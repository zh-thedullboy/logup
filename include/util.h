#pragma once

#include <cstdint>
#include <string>
#include <cassert>
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

// #include <iostream>
// int main(){
//     std::string a;
//     std::cin >> a;
//     size_t n = a.length();
//     std::cout<< convert_mask_to_u64(a, n) << '\n';
// }