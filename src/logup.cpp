#include "logup.h"
#include "goldilocks_quadratic_ext.h"
#include "mle_sumcheck.h"
#include "product_sumcheck.h"
#include "util.h"
#include <cassert>
#include <unordered_map>
#include <array>
#include <vector>
#include <random>

LogupProver::LogupProver(const table& f1, const table& f2, const table& t1, const table& t2):f1(f1), f2(f2), t1(t1), t2(t2) {
    calculate_multiplicities();
}

// 一定要用哈希表吗
void LogupProver::calculate_multiplicities(){
    size_t n = f1.size();
    size_t m = t1.size();
    assert(n == f2.size());
    assert(m == t2.size());

    // table c;
    std::vector<size_t> c1,c2;
    std::unordered_map<Goldilocks2::Element, size_t> freq_map1, freq_map2;
    for (size_t i = 0;i < n; ++i){
        freq_map1[f1[i]]++;
        freq_map2[f2[i]]++;
    }
    c1.reserve(t1.size());
    c2.reserve(t2.size());
    for(size_t i = 0;i < m; ++i){
        c1.push_back(freq_map1[t1[i]]);
        c2.push_back(freq_map2[t2[i]]);
    }




    
    assert(c1 == c2);
    c.reserve(c1.size());
    for (size_t i = 0;i < c1.size(); ++i){
        c.push_back(Goldilocks2::fromU64(c1[i]));
    }
    
    // print_table(t1);
    // print_table(t2);
    // print_table(f1);
    // print_table(f2);
    // print_table(c);
    // return MultilinearPolynomial(c);
}


void LogupProver::calculate_gh(const Goldilocks2::Element& gamma, const Goldilocks2::Element& lambda){
    h = c;
    g.resize(f1.size(), Goldilocks2::one());
    denomg.resize(g.size());
    denomh.resize(h.size());
    table inv(std::max(g.size(), h.size())); // store inverse of denomg and denomh

    for (size_t i = 0;i < f1.size(); ++i){
        Goldilocks2::Element tmp;
        Goldilocks2::mul(tmp, lambda, f2[i]);
        Goldilocks2::add(tmp, tmp, f1[i]);
        Goldilocks2::sub(denomg[i], gamma, tmp);
        // Goldilocks2::div(g[i], g[i], denomg[i]);
    }
    batch_inverse(inv, denomg);
    for(size_t i = 0;i < f1.size(); ++i){
        Goldilocks2::mul(g[i], g[i], inv[i]);
    }


    for (size_t i = 0;i < t1.size(); ++i){
        Goldilocks2::Element tmp;
        Goldilocks2::mul(tmp, lambda, t2[i]);
        Goldilocks2::add(tmp, tmp, t1[i]);
        Goldilocks2::sub(denomh[i], gamma, tmp);
        // Goldilocks2::div(h[i], h[i], denomh[i]);
    }
    batch_inverse(inv, denomh);
    for(size_t i = 0;i < t1.size(); ++i){
        Goldilocks2::mul(h[i], h[i], inv[i]);
    }


    // polyg = MultilinearPolynomial(g);
    polyg.emplace(g);
    polyh.emplace(h);
    // polyh = MultilinearPolynomial(h);
}

std::array<LogupProver::pcs, 4> LogupProver::commit_ft(){
    return {MultilinearPolynomial(f1), MultilinearPolynomial(f2), MultilinearPolynomial(t1), MultilinearPolynomial(t2)};
}

LogupProver::pcs LogupProver::commit_c(){
    return MultilinearPolynomial(c);
}

std::array<LogupProver::pcs, 2> LogupProver::commit_gh(){
    // return {MultilinearPolynomial(g), MultilinearPolynomial(h)};
    return {*polyg, *polyh};
}
std::array<LogupProver::pcs, 2> LogupProver::commit_denom(){
    // return {MultilinearPolynomial(g), MultilinearPolynomial(h)};
    return {MultilinearPolynomial(denomg), MultilinearPolynomial(denomh)};
}

std::array<sProver, 2> LogupProver::firstProvers(){
    return {sProver(*polyg), sProver(*polyh)};
}

std::array<pProver, 2> LogupProver::secondProvers(const std::vector<Goldilocks2::Element>& rg, const std::vector<Goldilocks2::Element>& rh){
    assert((1ull << rg.size()) == g.size());
    assert((1ull << rh.size()) == h.size());
    // assert(rh.size() == h.size());
    MultilinearPolynomial polyg(g);
    MultilinearPolynomial polyh(h);
    MultilinearPolynomial polydenomg(denomg);
    MultilinearPolynomial polydenomh(denomh);
    MultilinearPolynomial eqg = eq(polyg.get_num_vars(), rg);
    MultilinearPolynomial eqh = eq(polyh.get_num_vars(), rh);
    pProver prg(eqg, polyg, polydenomg);
    pProver prh(eqh, polyh, polydenomh);
    std::array<pProver, 2> provers = {prg, prh};
    return provers;
}


std::mt19937_64 LogupVerifier::gen(std::random_device{}());
std::uniform_int_distribution<uint64_t> LogupVerifier::dist(0, Goldilocks2::p - 1);


bool LogupVerifier::execute_logup(LogupProver& lpr){
    // inistialize the random generator
    // gen = std::mt19937_64(std::random_device{}());
    // constexpr uint64_t MODULUS = Goldilocks2::p;
    // dist = std::uniform_int_distribution<uint64_t>(1, MODULUS - 1);

    pcs c = lpr.commit_c();
    Goldilocks2::Element gamma = randnum();
    Goldilocks2::Element lambda = randnum();
    lpr.calculate_gh(gamma, lambda);

    std::array<pcs, 2> gh = lpr.commit_gh();
    // auto [g, h] = lpr.commit_gh();
    std::array<sProver, 2> firstProvers = lpr.firstProvers();
    Goldilocks2::Element sum = firstProvers[0].get_sum();
    assert(sum == firstProvers[1].get_sum());
    if(!sVerifier::execute_sumcheck(firstProvers[0], gh[0])){
        std::cout << "logup failed 0 \n";
        return false;
    }
    if(!sVerifier::execute_sumcheck(firstProvers[1], gh[1])){
        std::cout << "logup failed 1 \n";
        return false;
    }

    MultilinearPolynomial polyg = gh[0];
    MultilinearPolynomial polyh = gh[1];
    std::vector<Goldilocks2::Element> rg = randvec(gh[0].get_num_vars());
    std::vector<Goldilocks2::Element> rh = randvec(gh[1].get_num_vars());
    std::array<pProver, 2> secondProvers = lpr.secondProvers(rg, rh);
    assert(secondProvers[0].get_sum() == Goldilocks2::one());
    assert(secondProvers[1].get_sum() == c.evaluate(rh));

    // std::array<MultilinearPolynomial, 3> oracle
    MultilinearPolynomial eqg = eq(gh[0].get_num_vars(), rg);
    MultilinearPolynomial eqh = eq(gh[1].get_num_vars(), rh);
    std::array<pcs, 2> denom = lpr.commit_denom();

    if(!pVerifier::execute_sumcheck(secondProvers[0], {eqg, polyg, denom[0]})){
        std::cout << "logup failed 2 \n";
        return false;
    }
    if(!pVerifier::execute_sumcheck(secondProvers[1], {eqh, polyh, denom[1]})){
        std::cout << "logup failed 3 \n";
        return false;
    }
    // pVerifier::execute_sumcheck(secondProvers[0], {eqh, polyh, denom[1]});
    return true;
}
    
Goldilocks2::Element LogupVerifier::randnum(){
    return {Goldilocks::fromU64(dist(gen)), Goldilocks::fromU64(dist(gen))};
}

std::vector<Goldilocks2::Element> LogupVerifier::randvec(const uint64_t& n){
    std::vector<Goldilocks2::Element> rands;
    rands.reserve(n);
    for (uint64_t i = 0; i < n; ++i) {
        rands.push_back(randnum());
    }
    return rands;
}




// LogupProver::pcs LogupProver::commit_h() const{
//     return MultilinearPolynomial(h);
// }


// LogupProver::pcs LogupProver::commit_h(const Goldilocks2::Element& gamma, const Goldilocks2::Element& lambda) const{
//     table h = c;
//     for (size_t i = 0;i < t1.size(); ++i){
//         Goldilocks2::Element denom;
//         Goldilocks2::Element tmp;
//         Goldilocks2::mul(tmp, lambda, f2[i]);
//         Goldilocks2::add(tmp, tmp, f1[i]);
//         Goldilocks2::sub(denom, gamma, tmp);
//         Goldilocks2::div(g[i], g[i], denom);
//     }
//     return MultilinearPolynomial(g);
// }
