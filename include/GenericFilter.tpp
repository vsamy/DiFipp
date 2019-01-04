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

#include <limits>

namespace difi {

// Public static functions
template <typename T>
std::string GenericFilter<T>::filterStatus(FilterStatus status)
{
    switch (status) {
    case FilterStatus::NONE:
        return "Filter is uninitialized";
    case FilterStatus::READY:
        return "Filter is ready to be used";
    case FilterStatus::BAD_ORDER_SIZE:
        return "You try to initialize the filter with an order inferior or equal to 0 (window size for the moving average)";
    case FilterStatus::ALL_COEFF_MISSING:
        return "Filter has none of its coefficient initialized";
    case FilterStatus::A_COEFF_MISSING:
        return "Filter has its 'a' coefficients uninitialized";
    case FilterStatus::B_COEFF_MISSING:
        return "Filter has its 'b' coefficients uninitialized";
    case FilterStatus::BAD_FREQUENCY_VALUE:
        return "Filter has a received a frequency that is negative or equal to zero";
    case FilterStatus::BAD_CUTOFF_FREQUENCY:
        return "Filter has a received a bad cut-off frequency. It must be inferior to the sampling frequency";
    case FilterStatus::BAD_BAND_FREQUENCY:
        return "You try to initialize the filter with a bad combination of the frequency and bandwith, you must have fCenter > bw/2";
    default:
        return "I forgot to implement this error documentation";
    }
}

// Public functions

template <typename T>
T GenericFilter<T>::stepFilter(const T& data)
{
    assert(m_status == FilterStatus::READY);

    // Slide data (can't use SIMD, but should be small)
    for (Eigen::Index i = m_rawData.size() - 1; i > 0; --i)
        m_rawData(i) = m_rawData(i - 1);
    for (Eigen::Index i = m_filteredData.size() - 1; i > 0; --i)
        m_filteredData(i) = m_filteredData(i - 1);

    m_rawData[0] = data;
    m_filteredData[0] = 0;
    m_filteredData[0] = m_bCoeff.dot(m_rawData) - m_aCoeff.dot(m_filteredData);
    return m_filteredData[0];
}

template <typename T>
vectX_t<T> GenericFilter<T>::filter(const vectX_t<T>& data)
{
    vectX_t<T> results(data.size());
    if (!getFilterResults(results, data))
        return vectX_t<T>();

    return results;
}

template <typename T>
bool GenericFilter<T>::getFilterResults(Eigen::Ref<vectX_t<T>> results, const vectX_t<T>& data)
{
    assert(m_status == FilterStatus::READY);
    if (results.size() != data.size())
        return false;

    for (Eigen::Index i = 0; i < data.size(); ++i)
        results(i) = stepFilter(data(i));

    return true;
}

template <typename T>
void GenericFilter<T>::resetFilter()
{
    m_filteredData.setZero(m_aCoeff.size());
    m_rawData.setZero(m_bCoeff.size());
}

template <typename T>
template <typename T2>
void GenericFilter<T>::setCoeffs(T2&& aCoeff, T2&& bCoeff)
{
    static_assert(std::is_convertible_v<T2, vectX_t<T>>, "The coefficients types should be convertible to vectX_t<T>");

    if (!checkCoeffs(aCoeff, bCoeff))
        return;

    m_aCoeff = aCoeff;
    m_bCoeff = bCoeff;
    resetFilter();
    normalizeCoeffs();
}

template <typename T>
void GenericFilter<T>::getCoeffs(vectX_t<T>& aCoeff, vectX_t<T>& bCoeff) const
{
    aCoeff = m_aCoeff;
    bCoeff = m_bCoeff;
}

// Protected functions

template <typename T>
GenericFilter<T>::GenericFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff)
    : m_aCoeff(aCoeff)
    , m_bCoeff(bCoeff)
    , m_filteredData(aCoeff.size())
    , m_rawData(bCoeff.size())
{
    if (!checkCoeffs(aCoeff, bCoeff))
        return;

    resetFilter();
    normalizeCoeffs();
}

template <typename T>
void GenericFilter<T>::normalizeCoeffs()
{
    assert(m_status == FilterStatus::READY);

    T a0 = m_aCoeff(0);
    if (std::abs(a0 - T(1)) < std::numeric_limits<T>::epsilon())
        return;

    m_aCoeff /= a0;
    m_bCoeff /= a0;
}

template <typename T>
bool GenericFilter<T>::checkCoeffs(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff)
{
    m_status = FilterStatus::NONE;
    if (aCoeff.size() == 0)
        m_status = FilterStatus::A_COEFF_MISSING;
    else if (std::abs(aCoeff[0]) < std::numeric_limits<T>::epsilon())
        m_status = FilterStatus::BAD_A_COEFF;

    if (bCoeff.size() == 0)
        m_status = (m_status == FilterStatus::A_COEFF_MISSING ? FilterStatus::ALL_COEFF_MISSING : FilterStatus::B_COEFF_MISSING);

    if (m_status == FilterStatus::NONE)
        m_status = FilterStatus::READY;

    return m_status == FilterStatus::READY;
}

} // namespace difi