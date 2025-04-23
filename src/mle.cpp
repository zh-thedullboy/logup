#include "mle.h"
#include "util.h"
#include <string>

MultilinearPolynomial::MultilinearPolynomial(size_t num_vars):num_vars(num_vars){
    coeffs.resize(1ull << num_vars, Goldilocks2::zero());
}

void MultilinearPolynomial::set_term(std::string mask , const Goldilocks2::Element& c){
    coeffs[convert_mask_to_u64(mask, num_vars)] = c;
}

void MultilinearPolynomial::set_term(std::string mask , const uint64_t& c){
    uint64_t tmp[2] = {c, 0};
    coeffs[convert_mask_to_u64(mask, num_vars)] = Goldilocks2::fromU64(tmp);
}

Goldilocks2::Element MultilinearPolynomial::get_term(uint64_t mask) const{
    return coeffs[mask];
}

Goldilocks2::Element MultilinearPolynomial::evaluate(const std::vector<Goldilocks2::Element>& point) const{
    Goldilocks2::Element result = Goldilocks2::zero();
    for (uint64_t mask = 0; mask < coeffs.size(); ++mask) {
        Goldilocks2::Element term = coeffs[mask];
        for (size_t i = 0; i < num_vars; ++i) {
            if ((mask >> i) & 1) {
                // term *= point[i];
                Goldilocks2::mul(term, term, point[num_vars - i - 1]);
            }
        }
        // result += term;
        Goldilocks2::add(result, result, term);
    }
    return result;
}