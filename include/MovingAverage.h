// Copyright (c) 2019, Vincent SAMY
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the AIST nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include "DigitalFilter.h"
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
     */
    MovingAverage(int windowSize)
    {
        setWindowSize(windowSize);
    }
    /*! \brief Set the size of the moving average window. */
    void setWindowSize(int windowSize)
    {
        if (windowSize <= 0) {
            m_status = FilterStatus::BAD_ORDER_SIZE;
            return;
        }

        setCoeffs(vectX_t<T>::Constant(1, T(1)), vectX_t<T>::Constant(windowSize, T(1) / windowSize));
    }
    /*! \brief Get the size of the moving average window. */
    int windowSize() const noexcept { return bOrder(); }
};

} // namespace difi
