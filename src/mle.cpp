#include "mle.h"
#include "util.h"
#include <string>
#include <cassert>

MultilinearPolynomial::MultilinearPolynomial(size_t num_vars):num_vars(num_vars) {
    evaluations.resize(1ull << num_vars, Goldilocks2::zero());
}

MultilinearPolynomial::MultilinearPolynomial(const std::vector<Goldilocks2::Element>& evaluations):evaluations(evaluations) {
    size_t r = find_ceiling_log2(evaluations.size());
    num_vars = r;
}

MultilinearPolynomial::MultilinearPolynomial(const std::vector<uint64_t>& val_table){
    size_t r = find_ceiling_log2(val_table.size());
    num_vars = r;
    evaluations.resize(1ull << num_vars);
    for(size_t i = 0;i < val_table.size(); ++i){
        evaluations[i] = Goldilocks2::fromU64(val_table[i]);
    }
}

void MultilinearPolynomial::set_value(const std::string& mask, const Goldilocks2::Element& c){
    evaluations[convert_mask_to_u64(mask, num_vars)] = c;
}

void MultilinearPolynomial::set_value(const std::string& mask, const uint64_t& c){
    evaluations[convert_mask_to_u64(mask, num_vars)] = Goldilocks2::fromU64(c);
}

Goldilocks2::Element MultilinearPolynomial::eval_hypercube(uint64_t mask) const{
    return evaluations[mask];
}

Goldilocks2::Element MultilinearPolynomial::evaluate(const std::vector<Goldilocks2::Element>& point) const{
    // for denote purpose
    const std::vector<Goldilocks2::Element>& r = point;
    std::vector<Goldilocks2::Element> one_minus_r(num_vars);
    for (size_t i = 0; i < num_vars; ++i){
        Goldilocks2::sub(one_minus_r[i], Goldilocks2::one(), r[i]);
    }

    // construct lagrage bases
    std::vector<Goldilocks2::Element> lag_basis;
    lag_basis.resize(1ull << num_vars, Goldilocks2::one());
    // every round add a new bit to the highest digit, so the reversed ri should be utilized
    for(uint64_t i = 0;i < num_vars; ++i){
        for(uint64_t j = 0;j < (1ull << i); ++j){
            Goldilocks2::mul(lag_basis[j + (1ull << i)], lag_basis[j], r[num_vars - i - 1]);
            Goldilocks2::mul(lag_basis[j], lag_basis[j], one_minus_r[num_vars - i - 1]);
        }
    }

    Goldilocks2::Element result = Goldilocks2::zero();
    // directly uses lagrange bases rather than creates a new tmp var
    for (size_t i = 0; i < lag_basis.size(); ++i) {
        // Goldilocks2::Element tmp;
        Goldilocks2::mul(lag_basis[i], lag_basis[i], evaluations[i]);
        Goldilocks2::add(result, result, lag_basis[i]);
    }
    return result;
}

MultilinearPolynomial MultilinearPolynomial::operator+(const MultilinearPolynomial& g) const{
    assert(num_vars == g.get_num_vars());
    std::vector<Goldilocks2::Element> evs(evaluations.size()), evalg = g.get_eval_table();
    for(size_t i = 0;i < evaluations.size(); ++i){
        Goldilocks2::add(evs[i], evaluations[i], evalg[i]);
    }
    return MultilinearPolynomial(evs);
}


MultilinearPolynomial MultilinearPolynomial::operator-(const MultilinearPolynomial& g) const{
    assert(num_vars == g.get_num_vars());
    std::vector<Goldilocks2::Element> evs(evaluations.size()), evalg = g.get_eval_table();
    for(size_t i = 0;i < evaluations.size(); ++i){
        Goldilocks2::sub(evs[i], evaluations[i], evalg[i]);
    }
    return MultilinearPolynomial(evs);
}