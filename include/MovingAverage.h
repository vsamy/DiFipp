#pragma once

#include "GenericFilter.h"

namespace fratio {

template <typename T>
class MovingAverage : public GenericFilter<T> {
public:
    MovingAverage() = default;
    MovingAverage(size_t windowSize)
        : GenericFilter<T>({ 1 }, std::vector<T>(windowSize, T(1) / windowSize))
    {
    }

    void setWindowSize(size_t windowSize) { setCoeff({ 1 }, std::vector<T>(windowSize, 1 / windowSize)); }
    size_t windowSize() const noexcept { return m_bCoeff.size(); }
};

} // namespace fratio