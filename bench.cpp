#include "include/header"

/*
g++ -std=c++17 -O0 -o  bench bench.cpp ./src/* -lssl -lcrypto -lpthread -lgoldilocks -I./include -L./goldilocks/ -fopenmp -mavx2 -I./goldilocks/src -lgmp
*/

#include <vector>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <chrono>

#define iterations 1'000'000'000

void bench_operation(){
    Goldilocks::Element resb, base2, base1 = Goldilocks::fromU64(12898031213), base0 = Goldilocks::fromU64(314812375913);
    Goldilocks::Element base3 = Goldilocks::fromU64(11237412323), base4 = Goldilocks::fromU64(13951901234);
    
    
    Goldilocks2::Element rese, ext2, ext1 = {base0, base4}, ext0 = {base1, base3};

    auto start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < iterations; ++i){
        Goldilocks::mul(base1, base1, base0);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_ms = end - start;
    std::cout << "I took " << duration_ms.count() << " ms on the base field\n";
    std::cout << "dummy result: " << Goldilocks::toString(base1) << " \n";

    start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < iterations; ++i){
        Goldilocks2::mul(ext1, ext1, ext0);
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_ms_ext = end - start;
    std::cout << "I took " << duration_ms_ext.count() << " ms on the ext field\n";
    std::cout << "dummy result: " << Goldilocks2::toString(ext1) << " \n";
}

void bench_logup(const size_t& fsize){
    std::vector<Goldilocks2::Element> t1 = random_vec(1ull << 4);
    std::vector<Goldilocks2::Element> t2(t1.size());
    for(size_t i = 0;i < t1.size(); ++i){
        Goldilocks2::Element tmp = Goldilocks2::fromU64(2);
        Goldilocks2::mul(t2[i], t1[i], tmp);
    }

    std::vector<Goldilocks2::Element> f1(t1.size() << 3);
    std::vector<Goldilocks2::Element> f2(f1.size());
    srand(42);
    for(size_t i = 0;i < f1.size(); ++i){
        size_t r = rand() % t1.size();
        f1[i] = t1[r];
        f2[i] = t2[r];
    }

    LogupProver lpr(f1, f2, t1, t2);
    std::cout << LogupVerifier::execute_logup(lpr, 4, 64) << '\n';
}

int main(){
    bench_operation();
    
    
    return 0;
}