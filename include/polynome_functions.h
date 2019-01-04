#pragma once

#include "type_checks.h"
#include "typedefs.h"
#include <complex>

namespace fratio {

/*! \brief Compute polynome coefficients from roots.
 * 
 * This is done through Vieta's algorithm: \see https://en.wikipedia.org/wiki/Vieta%27s_formulas
 * \tparam T Floating type.
 */
template <typename T>
struct VietaAlgo {
    static_assert(std::is_arithmetic<internal::complex_sub_type_t<T>>::value, "This struct can only accept arithmetic types or complex.");

    /*! \brief Vieta's algorithm.
     * \note The function return the coefficients in the decreasing order: \f$a_n X^n + a_{n-1}X^{n-1} + ... + a1X + a0\f$.
     * \param roots Set of all roots of the polynome.
     * \return Coefficients of the polynome.
     */
    static vectX_t<T> polyCoeffFromRoot(const vectX_t<T>& roots);
};

template <typename T>
vectX_t<T> VietaAlgo<T>::polyCoeffFromRoot(const vectX_t<T>& roots)
{
    vectX_t<T> coeffs = vectX_t<T>::Zero(roots.size() + 1);
    coeffs(0) = T(1);
    for (Eigen::Index i = 0; i < roots.size(); ++i) {
        for (Eigen::Index k = i + 1; k > 0; --k) {
            coeffs(k) -= roots(i) * coeffs(k - 1);
        }
    }

    return coeffs;
}

} // namespace fratio
