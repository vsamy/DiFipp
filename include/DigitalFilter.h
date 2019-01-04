#pragma once

#include "GenericFilter.h"
#include "typedefs.h"

namespace fratio {

/*! \brief Basic digital filter.
 * 
 * This filter allows you to set any digital filter based on its coefficients.
 * \tparam T Floating type.
 */
template <typename T>
class DigitalFilter : public GenericFilter<T> {
public:
    /*! \brief Default uninitialized constructor. */
    DigitalFilter() = default;
    /*! \brief Constructor.
     * \param aCoeff Denominator coefficients of the filter in decreasing order.
     * \param bCoeff Numerator coefficients of the filter in decreasing order.
     */
    DigitalFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff)
        : GenericFilter<T>(aCoeff, bCoeff)
    {
    }
};

} // namespace fratio