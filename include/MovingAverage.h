#pragma once

#include "DigitalFilter.h"
#include "typedefs.h"

namespace fratio {

template <typename T>
class MovingAverage : public DigitalFilter<T> {
public:
    MovingAverage() = default;
    MovingAverage(int windowSize)
        : DigitalFilter<T>(vectX_t<T>::Constant(1, T(1)), vectX_t<T>::Constant(windowSize, T(1) / windowSize))
    {
    }

    void setWindowSize(int windowSize)
    {
        if (windowSize <= 0) {
            m_status = FilterStatus::BAD_ORDER_SIZE;
            return;
        }

        setCoeffs(vectX_t<T>::Constant(1, T(1)), vectX_t<T>::Constant(windowSize, T(1) / windowSize));
    }
    int windowSize() const noexcept { return bOrder(); }
};

} // namespace fratio
