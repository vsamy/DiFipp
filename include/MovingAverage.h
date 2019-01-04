#pragma once

#include "DigitalFilter.h"
#include "typedefs.h"

namespace fratio {

/*! \brief Moving average digital filter.
 * 
 * This is a specialization of a digital filter in order to use a moving average.
 * \tparam T Floating type.
 */
template <typename T>
class MovingAverage : public DigitalFilter<T> {
public:
    /*! \brief Default uninitialized constructor. */
    MovingAverage() = default;
    /*! \brief Constructor.
     * \param windowSize Size of the moving average window.
     */
    MovingAverage(int windowSize)
    {
        setWindowSize(windowSize);
    }
    /*! \brief Set the size of the moving average window. */
    void setWindowSize(int windowSize)
    {
        if (windowSize <= 0) {
            m_status = FilterStatus::BAD_ORDER_SIZE;
            return;
        }

        setCoeffs(vectX_t<T>::Constant(1, T(1)), vectX_t<T>::Constant(windowSize, T(1) / windowSize));
    }
    /*! \brief Get the size of the moving average window. */
    int windowSize() const noexcept { return bOrder(); }
};

} // namespace fratio
