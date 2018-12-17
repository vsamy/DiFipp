#pragma once

#include "GenericFilter.h"
#include "typedefs.h"

namespace fratio {

// https://www.mathworks.com/help/dsp/ug/how-is-moving-average-filter-different-from-an-fir-filter.html
template <typename T>
class DigitalFilter : public GenericFilter<T> {
public:
    DigitalFilter() = default;
    DigitalFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff)
        : GenericFilter<T>(aCoeff, bCoeff)
    {
    }
};

} // namespace fratio