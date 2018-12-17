#pragma once

#include "type_checks.h"
#include "typedefs.h"
#include <complex>

namespace fratio {

template <typename T>
struct VietaAlgo {
    static_assert(std::is_arithmetic<internal::complex_sub_type_t<T>>::value, "This struct can only accept arithmetic types or complex.");

    // Vieta's computation: https://en.wikipedia.org/wiki/Vieta%27s_formulas
    static vectX_t<T> polyCoeffFromRoot(const vectX_t<T>& poles);
};

template <typename T>
vectX_t<T> VietaAlgo<T>::polyCoeffFromRoot(const vectX_t<T>& poles)
{
    vectX_t<T> coeffs = vectX_t<T>::Zero(poles.size() + 1);
    coeffs(0) = T(1);
    for (Eigen::Index i = 0; i < poles.size(); ++i) {
        for (Eigen::Index k = i + 1; k > 0; --k) {
            coeffs(k) -= poles(i) * coeffs(k - 1);
        }
    }

    return coeffs;
}

} // namespace fratio
