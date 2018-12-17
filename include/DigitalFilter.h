#pragma once

#include "GenericFilter.h"
#include "typedefs.h"

namespace fratio {

template <typename T>
class DigitalFilter : public GenericFilter<T> {
public:
    DigitalFilter() = default;
    DigitalFilter(const Eigen::VectorX<T>& aCoeff, const Eigen::VectorX<T>& bCoeff)
        : GenericFilter<T>(aCoeff, bCoeff)
    {
    }
};

} // namespace fratio