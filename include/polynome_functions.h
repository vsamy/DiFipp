#pragma once

#include "type_checks.h"
#include <complex>
#include <vector>

namespace fratio {

template <typename T>
struct VietaAlgo {
    static_assert(std::is_arithmetic<internal::complex_sub_type_t<T>>::value, "This struct can only accept arithmetic types or complex.");

    // Vieta's computation: https://en.wikipedia.org/wiki/Vieta%27s_formulas
    static std::vector<T> polyCoeffFromRoot(const std::vector<T>& poles);
};

template <typename T>
std::vector<T> VietaAlgo<T>::polyCoeffFromRoot(const std::vector<T>& poles)
{
    std::vector<T> coeff(poles.size() + 1, 0);
    coeff[0] = 1;
    for (size_t i = 0; i < poles.size(); ++i) {
        for (size_t k = i + 1; k > 0; --k) {
            coeff[k] = coeff[k] - poles[i] * coeff[k - 1];
        }
    }

    return coeff;
}

} // namespace fratio