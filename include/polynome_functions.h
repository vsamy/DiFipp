// Copyright (c) 2019, Vincent SAMY
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies,
// either expressed or implied, of the FreeBSD Project.

#pragma once

#include "type_checks.h"
#include "typedefs.h"
#include <complex>

namespace difi {

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

} // namespace difi
