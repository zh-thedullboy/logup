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

LogupProver::LogupProver(const table_base& f_1, const table_base& f_2, const table_base& t_1, const table_base& t_2):f1(f_1), f2(f_2), t1(t_1), t2(t_2) {
    pad(f1, t1[0]);
    pad(f2, t2[0]);
    calculate_multiplicities();
}

// 一定要用哈希表吗
// void LogupProver::calculate_multiplicities(){
//     size_t n = f1.size();
//     size_t m = t1.size();
//     assert(n == f2.size());
//     assert(m == t2.size());

//     // table c;
//     std::vector<size_t> c1,c2;
//     std::unordered_map<Goldilocks::Element, size_t> freq_map1, freq_map2;
//     for (size_t i = 0;i < n; ++i){
//         freq_map1[f1[i]]++;
//         freq_map2[f2[i]]++;
//     }
//     c1.reserve(t1.size());
//     c2.reserve(t2.size());
//     for(size_t i = 0;i < m; ++i){
//         c1.push_back(freq_map1[t1[i]]);
//         c2.push_back(freq_map2[t2[i]]);
//     }
    
//     assert(c1 == c2);
//     c.reserve(c1.size());
//     for (size_t i = 0;i < c1.size(); ++i){
//         c.push_back(Goldilocks::fromU64(c1[i]));
//     }
    
//     // print_table(t1);
//     // print_table(t2);
//     // print_table(f1);
//     // print_table(f2);
//     // print_table(c);
//     // return MultilinearPolynomial(c);
// }


// void LogupProver::calculate_multiplicities(){
//     size_t n = f1.size();
//     size_t m = t1.size();
//     assert(n == f2.size());
//     assert(m == t2.size());

//     // table c;
//     std::vector<size_t> c1;
//     std::unordered_map<uint64_t, size_t> freq_map1, freq_map2;
//     for (size_t i = 0;i < n; ++i){
//         freq_map1[f1[i]]++;
//         freq_map2[f2[i]]++;
//     }
//     c1.reserve(t1.size());
//     c.reserve(t2.size());
//     for(size_t i = 0;i < m; ++i){
//         c1.push_back(freq_map1[t1[i]]);
//         c.push_back(freq_map2[t2[i]]);
//     }
    
//     assert(c1 == c);
// }


// in the real case, t1 and t2 are sorted vectors (ranges)
void LogupProver::calculate_multiplicities(){
    size_t n = f1.size();
    size_t m = t1.size();
    assert(n == f2.size());
    assert(m == t2.size());
    assert(is_power_of_2(m));
    
    c.resize(m, 0);

    bool flag = false;
    // don't calculate c twice.
    // instead, when you find f1[i] at t1[j], try to make sure f2[i] == t2[j]
    for(size_t i = 0; i < n; ++i){
        size_t idx = bisearch(t1, f1[i]);
        if(flag = flag || idx == t1.size() || f2[i] != t2[idx]) break;
        c[idx] += 1;
    }
    assert(!flag);
}


void LogupProver::calculate_gh(const Goldilocks2::Element& gamma, const Goldilocks2::Element& lambda){
    // h = c;

    h.resize(c.size());
    for(size_t i = 0;i < c.size(); ++i){
        h[i] = {c[i], Goldilocks::zero()};
    }
    g.resize(f1.size(), Goldilocks2::one());
    denomg.resize(g.size());
    denomh.resize(h.size());
    table_ext inv(std::max(g.size(), h.size())); // store inverse of denomg and denomh

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


    polyg.emplace(g);
    polyh.emplace(h);
}

std::array<LogupDef::pcs_base, 4> LogupProver::commit_ft(const uint64_t& rho_inv){
    ligeroProver_base 
                prf1(f1, rho_inv),
                prf2(f2, rho_inv),
                prt1(t1, rho_inv),
                prt2(t2, rho_inv);
    return {prf1.commit(), prf2.commit(), prt1.commit(), prt2.commit()};
}

LogupDef::pcs_base LogupProver::commit_c(const uint64_t& rho_inv){
    // return MultilinearPolynomial(c);
    ligeroProver_base pr(c, rho_inv);
    return pr.commit();
}

std::array<LogupDef::pcs_ext, 2> LogupProver::commit_gh(const uint64_t& rho_inv){
    ligeroProver_ext prg(g, rho_inv), prh(h, rho_inv);
    return {prg.commit(), prh.commit()};
}


std::array<sProver, 2> LogupProver::firstProvers(){
    return {sProver(*polyg), sProver(*polyh)};
}

std::array<pProver, 2> LogupProver::secondProvers(const std::vector<Goldilocks2::Element>& rg, const std::vector<Goldilocks2::Element>& rh){
    assert((1ull << rg.size()) == g.size());
    assert((1ull << rh.size()) == h.size());
    // assert(rh.size() == h.size());
    MultilinearPolynomial polydenomg(denomg);
    MultilinearPolynomial polydenomh(denomh);
    MultilinearPolynomial eqg = eq((*polyg).get_num_vars(), rg);
    MultilinearPolynomial eqh = eq((*polyh).get_num_vars(), rh);
    pProver prg(eqg, *polyg, polydenomg);
    pProver prh(eqh, *polyh, polydenomh);
    std::array<pProver, 2> provers = {prg, prh};
    return provers;
}


std::mt19937_64 LogupVerifier::gen(std::random_device{}());
std::uniform_int_distribution<uint64_t> LogupVerifier::dist(0, Goldilocks2::p - 1);

bool LogupVerifier::execute_logup(LogupProver& lpr, const uint64_t& rho_inv, const size_t& sec_param){
    auto ft = lpr.commit_ft(rho_inv);
    for(auto pc: ft){
        if(!ligeroVerifier::check_commit(pc, sec_param)) return false;
    }
    alert("f,t commited");

    auto c = lpr.commit_c(rho_inv);
    if(!ligeroVerifier::check_commit(c, sec_param)) return false;
    alert("c commited");

    Goldilocks2::Element gamma = randnum();
    Goldilocks2::Element lambda = randnum();
    lpr.calculate_gh(gamma, lambda);
    auto gh = lpr.commit_gh(rho_inv);
    for(auto pc: gh){
        if(!ligeroVerifier::check_commit(pc, sec_param)) return false;
    }
    alert("g,h commited");

    auto pcsg = gh[0], pcsh = gh[1];

    std::array<sProver, 2> firstProvers = lpr.firstProvers();
    Goldilocks2::Element sum = firstProvers[0].get_sum();
    assert(sum == firstProvers[1].get_sum());
    if(!sVerifier::execute_sumcheck(firstProvers[0], pcsg, sec_param)){
        std::cout << "logup failed 0 \n";
        return false;
    }
    alert("sumcheck 1 / 4 finished");
    if(!sVerifier::execute_sumcheck(firstProvers[1], pcsh, sec_param)){
        std::cout << "logup failed 1 \n";
        return false;
    }
    alert("sumcheck 2 / 4 finished");

    const size_t numvar_g = find_ceiling_log2(pcsg.num_cols * pcsg.num_rows);
    const size_t numvar_h = find_ceiling_log2(pcsh.num_cols * pcsh.num_rows);
    std::vector<Goldilocks2::Element> rg = randvec(numvar_g);
    std::vector<Goldilocks2::Element> rh = randvec(numvar_h);

    std::array<pProver, 2> secondProvers = lpr.secondProvers(rg, rh);
    assert(secondProvers[0].get_sum() == Goldilocks2::one());
    assert(secondProvers[1].get_sum() == ligeroVerifier::open(c, rh, sec_param));


    MultilinearPolynomial eqg = eq(numvar_g, rg);
    MultilinearPolynomial eqh = eq(numvar_h, rh);
    auto pcsf1 = ft[0], pcsf2 = ft[1], pcst1 = ft[2], pcst2 = ft[3];
    if(!pVerifier::execute_logup_sumcheck(secondProvers[0], eqg, pcsg, pcsf1, pcsf2, gamma, lambda, sec_param)){
        std::cout << "logup failed 2 \n";
        return false;
    }
    alert("sumcheck 3 / 4 finished");

    if(!pVerifier::execute_logup_sumcheck(secondProvers[1], eqh, pcsh, pcst1, pcst2, gamma, lambda, sec_param)){
        std::cout << "logup failed 3 \n";
        return false;
    }
    alert("sumcheck 4 / 4 finished");
    alert("logup finished");
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




// LogupDef::pcs LogupProver::commit_h() const{
//     return MultilinearPolynomial(h);
// }


// LogupDef::pcs LogupProver::commit_h(const Goldilocks2::Element& gamma, const Goldilocks2::Element& lambda) const{
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
