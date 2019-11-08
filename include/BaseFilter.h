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

#include "gsl/gsl_assert.h"
#include "type_checks.h"
#include "typedefs.h"
#include <string>

namespace difi {

// TODO: noexcept(Function of gsl variable)
// TODO: constructor with universal refs
// TODO: setCoeffs with universal refs

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
template <typename T, typename Derived>
class BaseFilter {
    static_assert(std::is_floating_point<T>::value && !std::is_const<T>::value, "Only accept non-complex floating point types.");
    friend Derived;

public:
    /*! \brief Reset the data and filtered data. */
    void resetFilter() noexcept { derived().resetFilter(); };

    /*!< \brief Return the filter type */
    FilterType type() const noexcept { return m_type; }
    /*! \brief Get how far back is the filtered value.
     * The filtered value is computed at time \f$tc=t - center()*\Delta T\f$ where \f$t\f$ is the current time
     * \f$\Delta T\f$ is your sampling time. 
     * \note For FilterType::Backward filter, the function returns 0.
     */
    Eigen::Index center() const noexcept { return (m_type == FilterType::Backward ? 0 : ((m_bCoeff.size() - 1) / 2)); }
    /*! \brief Get digital filter coefficients.
     * 
     * It will automatically resize the given vectors.
     * \param[out] aCoeff Denominator coefficients of the filter in decreasing order.
     * \param[out] bCoeff Numerator coefficients of the filter in decreasing order.
     */
    void getCoeffs(vectX_t<T>& aCoeff, vectX_t<T>& bCoeff) const noexcept;
    /*! \brief Return coefficients of the denominator polynome. */
    const vectX_t<T>& aCoeff() const noexcept { return m_aCoeff; }
    /*! \brief Return coefficients of the numerator polynome. */
    const vectX_t<T>& bCoeff() const noexcept { return m_bCoeff; }
    /*! \brief Return the order the denominator polynome order of the filter. */
    Eigen::Index aOrder() const noexcept { return m_aCoeff.size(); }
    /*! \brief Return the order the numerator polynome order of the filter. */
    Eigen::Index bOrder() const noexcept { return m_bCoeff.size(); }
    /*! \brief Return the initialization state of the filter0 */
    bool isInitialized() const noexcept { return m_isInitialized; }
    /*! \brief Set type of filter (one-sided or centered)
     * 
     * \param type The filter type.
     * \warning bCoeff must be set before.
     */
    void setType(FilterType type);
    /*! \brief Set the new coefficients of the filters.
     *
     * It awaits a universal reference.
     * \param aCoeff Denominator coefficients of the filter in decreasing order.
     * \param bCoeff Numerator coefficients of the filter in decreasing order.
     */
    void setCoeffs(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff);

protected:
    /*! \brief Normalized the filter coefficients such that aCoeff(0) = 1. */
    void normalizeCoeffs();
    /*! \brief Check for bad coefficients.
     * 
     * Set the filter status to ready is everything is fine.
     * \param aCoeff Denominator coefficients of the filter.
     * \param bCoeff Numerator coefficients of the filter.
     * \return True if the filter status is set on READY.
     */
    bool checkCoeffs(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff);

private:
    /*! \brief Default uninitialized constructor. */
    BaseFilter() = default;
    /*! \brief Constructor.
     * \param aCoeff Denominator coefficients of the filter in decreasing order.
     * \param bCoeff Numerator coefficients of the filter in decreasing order.
     * \param type Type of the filter.
     */
    BaseFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff, FilterType type = FilterType::Backward);
    /*! \brief Default destructor. */
    virtual ~BaseFilter() = default;

    Derived& derived() noexcept { return *static_cast<Derived*>(this); }
    const Derived& derived() const noexcept { return *static_cast<const Derived*>(this); }

private:
    FilterType m_type; /*!< Type of filter. Default is FilterType::Backward. */
    bool m_isInitialized = false; /*!< Initialization state of the filter. Default is false */
    vectX_t<T> m_aCoeff; /*!< Denominator coefficients of the filter */
    vectX_t<T> m_bCoeff; /*!< Numerator coefficients of the filter */
    vectX_t<T> m_filteredData; /*!< Last set of filtered data */
    vectX_t<T> m_rawData; /*!< Last set of non-filtered data */
};

} // namespace difi

#include "BaseFilter.tpp"
