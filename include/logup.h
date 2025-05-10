#pragma once

#include "goldilocks_quadratic_ext.h"
#include "mle_sumcheck.h"
#include "product_sumcheck.h"
#include "mle.h"
#include "ligero.h"
#include <vector>
#include <array>
#include <optional>
#include <random>

namespace LogupDef{
    typedef ligeropcs pcs;
}

class LogupProver{
public:
    using table = std::vector<Goldilocks2::Element>;
    // should be replaced with a pcs
private:
    table f1, f2, t1, t2, c, g, h;
    std::optional<MultilinearPolynomial> polyg, polyh;
    // intermediate tables used for last 2 sumchecks
    table denomg, denomh;
    // std::optional<sProver> sp;
    // std::optional<pProver> pp;
public:
    LogupProver(const table& f1, const table& f2, const table& t1, const table& t2);
    void calculate_multiplicities();
    void calculate_gh(const Goldilocks2::Element& gamma, const Goldilocks2::Element& lambda);
    LogupDef::pcs commit_c(const uint64_t& rho_inv);
    std::array<LogupDef::pcs, 4> commit_ft(const uint64_t& rho_inv);
    std::array<LogupDef::pcs, 2> commit_gh(const uint64_t& rho_inv);
    std::array<sProver, 2> firstProvers();
    std::array<pProver, 2> secondProvers(const std::vector<Goldilocks2::Element>& rg, const std::vector<Goldilocks2::Element>& rh);
};

class LogupVerifier{
public:
    static bool execute_logup(LogupProver& lpr, const uint64_t& rho_inv, const size_t& sec_param);
private:
    static std::mt19937_64 gen;
    static std::uniform_int_distribution<uint64_t> dist;
    static Goldilocks2::Element randnum();
    static std::vector<Goldilocks2::Element> randvec(const uint64_t& n);
};