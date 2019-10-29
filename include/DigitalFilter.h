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

#include "GenericFilter.h"
#include "typedefs.h"

namespace difi {

/*! \brief Basic digital filter.
 * 
 * This filter allows you to set any digital filter based on its coefficients.
 * \tparam T Floating type.
 */
template <typename T>
class DigitalFilter : public GenericFilter<T> {
public:
    /*! \brief Default uninitialized constructor. */
    DigitalFilter() = default;
    /*! \brief Constructor.
     * \param aCoeff Denominator coefficients of the filter in decreasing order.
     * \param bCoeff Numerator coefficients of the filter in decreasing order.
     */
    DigitalFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff)
        : GenericFilter<T>(aCoeff, bCoeff)
    {
    }

    void setCoefficients(vectX_t<T>&& aCoeff, vectX_t<T>&& bCoeff)
    {
        setCoeffs(std::forward(aCoeff), std::forward(bCoeff));
    }
};

/*! \brief Basic centered digital filter.
 * 
 * This filter allows you to set any centered digital filter based on its coefficients.
 * \tparam T Floating type.
 */
template <typename T>
class CenteredDigitalFilter : public GenericFilter<T> {
public:
    /*! \brief Default uninitialized constructor. */
    CenteredDigitalFilter() = default;
    /*! \brief Constructor.
     * \param aCoeff Denominator coefficients of the filter in decreasing order.
     * \param bCoeff Numerator coefficients of the filter in decreasing order.
     */
    CenteredDigitalFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff)
        : GenericFilter<T>(aCoeff, bCoeff, Type::Centered)
    {
    }
    void setCoefficients(vectX_t<T>&& aCoeff, vectX_t<T>&& bCoeff)
    {
        setCoeffs(std::forward(aCoeff), std::forward(bCoeff));
        setType(Type::Centered);
    }
};

} // namespace difi