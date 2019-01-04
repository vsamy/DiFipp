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

#include <Eigen/Core>

namespace difi {

template <typename T>
using vectX_t = Eigen::Matrix<T, Eigen::Dynamic, 1>; /*!< Eigen column-vector */

template <typename T>
using vectXc_t = vectX_t<std::complex<T>>; /*!< Eigen complex column-vector */

/*! \brief Filter status */
enum class FilterStatus {
    // Generic filter
    NONE, /*!< Filter has not yet been initialized */
    READY, /*!< Filter is ready to process data */
    BAD_ORDER_SIZE, /*!< Order of the filter is bad */
    BAD_A_COEFF, /*!< Denominator coefficients of the filter is bad */
    A_COEFF_MISSING, /*!< Denominator coefficients have not been set */
    B_COEFF_MISSING, /*!< Numerator coefficients have not been set */
    ALL_COEFF_MISSING = A_COEFF_MISSING | B_COEFF_MISSING, /*!< Coefficients have not been set */

    // Butterworth filter
    BAD_FREQUENCY_VALUE, /*!< Given frequency is bad */
    BAD_CUTOFF_FREQUENCY, /*!< Given ctu-off frequency is bad */
    BAD_BAND_FREQUENCY /*!< Given band frequency is bad */
};

} // namespace difi