#pragma once

#include <Eigen/Core>

namespace fratio {

template <typename T>
using vectX_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;

enum class FilterStatus {
    // Generic filter
    NONE,
    READY,
    BAD_ORDER_SIZE,
    BAD_A_COEFF,
    A_COEFF_MISSING,
    B_COEFF_MISSING,
    ALL_COEFF_MISSING = A_COEFF_MISSING | B_COEFF_MISSING,

    // Butterworth filter
    BAD_FREQUENCY_VALUE,
    BAD_CUTOFF_FREQUENCY
};

} // namespace fratio