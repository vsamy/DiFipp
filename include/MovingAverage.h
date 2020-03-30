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

#include "DigitalFilter.h"
#include "gsl/gsl_assert.h"
#include "typedefs.h"

namespace difi {

/*! \brief Moving average digital filter.
 * 
 * This is a specialization of a digital filter in order to use a moving average.
 * \tparam T Floating type.
 */
template <typename T>
class MovingAverage : public DigitalFilter<T> {
public:
    /*! \brief Default uninitialized constructor. */
    MovingAverage() = default;
    /*! \brief Constructor.
     * \param windowSize Size of the moving average window.
     * \param type Type of the filter.
     */
    MovingAverage(int windowSize, FilterType type = FilterType::Backward)
    {
        setWindowSize(windowSize);
        this->setType(type);
    }
    /*! \brief Set the size of the moving average window. */
    void setWindowSize(int windowSize)
    {
        Expects(windowSize > 0);
        this->setCoeffs(vectX_t<T>::Constant(1, T(1)), vectX_t<T>::Constant(windowSize, T(1) / static_cast<T>(windowSize)));
    }
    /*! \brief Get the size of the moving average window. */
    int windowSize() const noexcept { return this->bOrder(); }
};

} // namespace difi
