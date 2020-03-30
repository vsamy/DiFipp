// Copyright (c) 2019, Vincent SAMY
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies,
// either expressed or implied, of the FreeBSD Project.

#pragma once

#include "BaseFilter.h"

namespace difi {

template <typename T>
class GenericFilter : public BaseFilter<T, GenericFilter<T>> {
    using Base = BaseFilter<T, GenericFilter<T>>;
    using Base::m_isInitialized;
    using Base::m_aCoeff;
    using Base::m_bCoeff;
    using Base::m_rawData;
    using Base::m_filteredData;

public:
    /*! \brief Filter a new data.
     * 
     * This function is practical for online application that does not know the whole signal in advance.
     * \param data New data to filter.
     * \return Filtered data.
     */
    T stepFilter(const T& data);
    /*! \brief Filter a signal.
     * 
     * Filter all data given by the signal.
     * \param data Signal.
     * \return Filtered signal.
     */
    vectX_t<T> filter(const vectX_t<T>& data);

    void resetFilter() noexcept;

protected:
    GenericFilter() = default;
    GenericFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff, FilterType type = FilterType::Backward)
        : Base(aCoeff, bCoeff, type)
    {}
};

template <typename T>
class TVGenericFilter : public BaseFilter<T, TVGenericFilter<T>> {
    using Base = BaseFilter<T, TVGenericFilter<T>>;
    using Base::m_isInitialized;
    using Base::m_aCoeff;
    using Base::m_bCoeff;
    using Base::m_rawData;
    using Base::m_filteredData;

public:
    /*! \brief Filter a new data.
     * 
     * This function is practical for online application that does not know the whole signal in advance.
     * \param data New data to filter.
     * \return Filtered data.
     */
    T stepFilter(const T& time, const T& data);
    /*! \brief Filter a signal.
     * 
     * Filter all data given by the signal.
     * \param data Signal.
     * \return Filtered signal.
     */
    vectX_t<T> filter(const vectX_t<T>& data, const vectX_t<T>& time);

    void resetFilter() noexcept;

protected:
    TVGenericFilter() = default;
    TVGenericFilter(size_t differentialOrder, const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff, FilterType type = FilterType::Backward)
        : Base()
        , m_diffOrder(differentialOrder)
    {
        Expects(differentialOrder >= 1);
        this->setCoeffs(aCoeff, bCoeff);
        this->setType(type);
    }

private:
    size_t m_diffOrder = 1;
    vectX_t<T> m_timers;
};

} // namespace difi

#include "GenericFilter.tpp"