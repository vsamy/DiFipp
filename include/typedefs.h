#pragma once

#include <Eigen/Core>

namespace fratio {

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

} // namespace fratio