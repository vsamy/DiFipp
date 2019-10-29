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
class GenericFilter : BaseFilter<T, GenericFilter<T>> {
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

    void resetFilter() noexcept 

protected:
    GenericFilter() = default;
    GenericFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff, Type type = Type::Forward)
        : BaseFilter(aCoeff, bCoeff, type)
    {}
};

template <typename T>
class TVGenericFilter : BaseFilter<T, GenericFilter<T>> {
public:
    /*! \brief Filter a new data.
     * 
     * This function is practical for online application that does not know the whole signal in advance.
     * \param data New data to filter.
     * \return Filtered data.
     */
    T stepFilter(const T& data, const T& time);
    /*! \brief Filter a signal.
     * 
     * Filter all data given by the signal.
     * \param data Signal.
     * \return Filtered signal.
     */
    vectX_t<T> filter(const vectX_t<T>& data, const vectX_t<T>& time);

protected:
    GenericFilter() = default;
    GenericFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff, Type type = Type::Forward)
        : BaseFilter(aCoeff, bCoeff, type)
    {
        m_timers.resize(m_bCoeffs.size());
        m_timeDiffs.resize(m_bCoeffs.size());
    }

private:
    vectX_t<T> m_timers;
    vectX_t<T> m_timeDiffs;
};

} // namespace difi

#include "GenericFilter.tpp"