#pragma once

#include "type_checks.h"
#include "typedefs.h"
#include <complex>

namespace fratio {

template <typename T>
struct VietaAlgo {
    static_assert(std::is_arithmetic<internal::complex_sub_type_t<T>>::value, "This struct can only accept arithmetic types or complex.");

    // Vieta's computation: https://en.wikipedia.org/wiki/Vieta%27s_formulas
    static Eigen::VectorX<T> polyCoeffFromRoot(const Eigen::VectorX<T>& poles);
};

template <typename T>
Eigen::VectorX<T> VietaAlgo<T>::polyCoeffFromRoot(const Eigen::VectorX<T>& poles)
{
    Eigen::VectorX<T> coeffs = Eigen::VectorX<T>::Zero(poles.size() + 1);
    coeffs(0) = T(1);
    for (size_t i = 0; i < poles.size(); ++i) {
        for (size_t k = i + 1; k > 0; --k) {
            coeffs(k) -= poles(i) * coeffs(k - 1);
        }
    }
    // Check for equation c(k) = sum(i=k-1, poles.size() : p(i)) * c(k-1), k>=1
    // size_t pSize = poles.size();
    // for (size_t k = 1; k < coeffs.size(); ++k)
    //     coeffs(k) -= poles.tail(pSize - (k - 1)).sum() * coeffs(k - 1);

    // Maybe better
    // T sum = poles.sum();
    // for (size_t k = 1; k < coeffs.size(); ++k) {
    //     coeffs(k) -= sum * coeffs(k - 1);
    //     sum -= poles(k - 1);
    // } 

    return coeffs;
}

} // namespace fratio
