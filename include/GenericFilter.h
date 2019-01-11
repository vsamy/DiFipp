// Copyright (c) 2019, Vincent SAMY
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the AIST nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include "type_checks.h"
#include "typedefs.h"
#include <stddef.h>
#include <string>

namespace difi {

/*! \brief Low-level filter.
 * 
 * It creates the basic and common functions of all linear filter that can written as a digital filter.
 * This class can not be instantiated directly.
 * 
 * \warning In Debug mode, all functions may throw if a filter is badly initialized.
 * This not the case in Realese mode.
 * 
 * \tparam T Floating type.
 */
template <typename T>
class GenericFilter {
    static_assert(std::is_floating_point<T>::value && !std::is_const<T>::value, "Only accept non-complex floating point types.");

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
    /*! \brief Filter a signal and store in a user-defined Eigen vector.
     * 
     * Useful if the length of the signal is known in advance.
     * \param[out] results Filtered signal.
     * \param data Signal.
     * \return False if vector's lengths do not match.
     */
    void getFilterResults(Eigen::Ref<vectX_t<T>> results, const vectX_t<T>& data);
    /*! \brief Reset the data and filtered data. */
    void resetFilter() noexcept;
    /*! \brief Get digital filter coefficients.
     * 
     * It will automatically resize the given vectors.
     * \param[out] aCoeff Denominator coefficients of the filter in decreasing order.
     * \param[out] bCoeff Numerator coefficients of the filter in decreasing order.
     */
    void getCoeffs(vectX_t<T>& aCoeff, vectX_t<T>& bCoeff) const noexcept;
    /*! \brief Return the order the denominator polynome order of the filter. */
    Eigen::Index aOrder() const noexcept { return m_aCoeff.size(); }
    /*! \brief Return the order the numerator polynome order of the filter. */
    Eigen::Index bOrder() const noexcept { return m_bCoeff.size(); }
    /*! \brief Return the initialization state of the filter0 */
    bool isInitialized() const noexcept { return m_isInitialized; }

protected:
    /*! \brief Default uninitialized constructor. */
    GenericFilter() = default;
    /*! \brief Constructor.
     * \param aCoeff Denominator coefficients of the filter in decreasing order.
     * \param bCoeff Numerator coefficients of the filter in decreasing order.
     */
    GenericFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff);
    /*! \brief Default destructor. */
    virtual ~GenericFilter() = default;

    /*! \brief Set the new coefficients of the filters.
     *
     * It awaits a universal reference.
     * \param aCoeff Denominator coefficients of the filter in decreasing order.
     * \param bCoeff Numerator coefficients of the filter in decreasing order.
     */
    template <typename T2>
    void setCoeffs(T2&& aCoeff, T2&& bCoeff);
    /*! \brief Normalized the filter coefficients such that aCoeff(0) = 1. */
    void normalizeCoeffs();
    /*! \brief Check for bad coefficients.
     * 
     * Set the filter status to ready is everything is fine.
     * \param aCoeff Denominator coefficients of the filter.
     * \param bCoeff Numerator coefficients of the filter.
     * \return True if the filter status is set on READY.
     */
    void checkCoeffs(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff);

private:
    bool m_isInitialized = false; /*!< Initialization state of the filter. Default is false */
    vectX_t<T> m_aCoeff; /*!< Denominator coefficients of the filter */
    vectX_t<T> m_bCoeff; /*!< Numerator coefficients of the filter */
    vectX_t<T> m_filteredData; /*!< Last set of filtered data */
    vectX_t<T> m_rawData; /*!< Last set of non-filtered data */
};

} // namespace difi

#include "GenericFilter.tpp"
