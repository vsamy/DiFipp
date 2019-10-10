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

#include "gsl/gsl_assert.h"
#include "type_checks.h"
#include "typedefs.h"
#include <limits>

namespace difi {

/*! \brief Transform an analog signal to a discrete signal and vice versa.
 * 
 * \see https://en.wikipedia.org/wiki/Bilinear_transform
 * \tparam T Floating (complex) types.
 */
template <typename T>
struct BilinearTransform {
    using SubType = internal::complex_sub_type_t<T>; /*!< Sub-type of the complex if T is complex, T otherwise */
    static_assert(std::is_floating_point<SubType>::value, "This struct can only accept floating point types (real and complex).");

    /*! \brief Transformation from analog to discrete.
     * \param fs Sampling frequency.
     * \param sPlanePole Analog data.
     * \param[out] zPlanePole Resulting discrete data.
     */
    static void SToZ(SubType fs, const T& sPlanePole, T& zPlanePole);
    /*! \brief Transformation from analog to discrete.
     * \param fs Sampling frequency.
     * \param sPlanePole Analog signal.
     * \param[out] zPlanePole Resulting discrete signal.
     */
    static void SToZ(SubType fs, const vectX_t<T>& sPlanePoles, Eigen::Ref<vectX_t<T>>& zPlanePoles); // Can be optimized maybe
    /*! \brief Transformation from discrete to analog.
     * \param fs Sampling frequency.
     * \param zPlanePole Discrete data.
     * \param[out] sPlanePole Resulting analog data.
     */
    static void ZToS(SubType fs, const T& zPlanePole, T& sPlanePole);
    /*! \brief Transformation from discrete to analog.
     * \param fs Sampling frequency.
     * \param zPlanePole Discrete signal.
     * \param[out] sPlanePole Resulting analog signal.
     */
    static void ZToS(SubType fs, const vectX_t<T>& zPlanePoles, Eigen::Ref<vectX_t<T>>& sPlanePoles); // Can be optimized maybe
};

template <typename T>
void BilinearTransform<T>::SToZ(SubType fs, const T& sPlanePole, T& zPlanePole)
{
    Expects(std::abs(2 * fs - sPlanePole) > std::numeric_limits<SubType>::epsilon()); // Divide-by-zero otherwise
    T scalePole = sPlanePole / (2 * fs);
    zPlanePole = (T(1) + scalePole) / (T(1) - scalePole);
}

template <typename T>
void BilinearTransform<T>::SToZ(SubType fs, const vectX_t<T>& sPlanePoles, Eigen::Ref<vectX_t<T>>& zPlanePoles)
{
    Expects(sPlanePoles.size() == zPlanePoles.size());
    for (Eigen::Index k = 0; k < sPlanePoles.size(); ++k)
        SToZ(fs, sPlanePoles(k), zPlanePoles(k));
}

template <typename T>
void BilinearTransform<T>::ZToS(SubType fs, const T& zPlanePole, T& sPlanePole)
{
    Expects(std::abs(T(1) + zPlanePole) > std::numeric_limits<SubType>::epsilon()); // Divide-by-zero otherwise
    T invPole = T(1) / zPlanePole;
    sPlanePole = 2 * fs * (T(1) - invPole) / (T(1) + invPole);
}

template <typename T>
void BilinearTransform<T>::ZToS(SubType fs, const vectX_t<T>& zPlanePoles, Eigen::Ref<vectX_t<T>>& sPlanePoles)
{
    Expects(sPlanePoles.size() == zPlanePoles.size());
    for (Eigen::Index k = 0; k < sPlanePoles.size(); ++k)
        ZToS(fs, zPlanePoles(k), sPlanePoles(k));
}

} // namespace difi