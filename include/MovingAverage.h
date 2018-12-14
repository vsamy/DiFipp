#pragma once

#include "DigitalFilter.h"
#include "typedefs.h"

namespace fratio {

template <typename T>
class MovingAverage : public DigitalFilter<T> {
public:
    MovingAverage() = default;
    MovingAverage(size_t windowSize)
        : DigitalFilter<T>(Eigen::VectorX<T>::Constant(1, T(1)), Eigen::VectorX<T>::Constant(windowSize, T(1) / windowSize))
    {
    }

    void setWindowSize(size_t windowSize) { setCoeffs(Eigen::VectorX<T>::Constant(1, T(1)), Eigen::VectorX<T>::Constant(windowSize, T(1) / windowSize)); }
    size_t windowSize() const noexcept { return bOrder(); }
};

} // namespace fratio
