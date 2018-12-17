#pragma once

#include <Eigen/Core>

namespace Eigen {

template <typename T>
using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

} // namespace Eigen

namespace fratio {

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