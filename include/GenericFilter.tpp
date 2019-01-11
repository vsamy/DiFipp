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

// Public functions

template <typename T>
T GenericFilter<T>::stepFilter(const T& data)
{
    Expects(m_isInitialized);

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
    getFilterResults(results, data);
    return results;
}

template <typename T>
void GenericFilter<T>::getFilterResults(Eigen::Ref<vectX_t<T>> results, const vectX_t<T>& data)
{
    Expects(m_isInitialized);
    Expects(results.size() == data.size());
    for (Eigen::Index i = 0; i < data.size(); ++i)
        results(i) = stepFilter(data(i));
}

template <typename T>
void GenericFilter<T>::resetFilter() noexcept
{
    m_filteredData.setZero(m_aCoeff.size());
    m_rawData.setZero(m_bCoeff.size());
}

template <typename T>
template <typename T2>
void GenericFilter<T>::setCoeffs(T2&& aCoeff, T2&& bCoeff)
{
    static_assert(std::is_convertible_v<T2, vectX_t<T>>, "The coefficients types should be convertible to vectX_t<T>");

    checkCoeffs(aCoeff, bCoeff);
    m_aCoeff = aCoeff;
    m_bCoeff = bCoeff;
    normalizeCoeffs();
    resetFilter();
    m_isInitialized = true;
}

template <typename T>
void GenericFilter<T>::getCoeffs(vectX_t<T>& aCoeff, vectX_t<T>& bCoeff) const noexcept
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
    checkCoeffs(aCoeff, bCoeff);
    normalizeCoeffs();
    resetFilter();
    m_isInitialized = true;
}

template <typename T>
void GenericFilter<T>::normalizeCoeffs()
{
    T a0 = m_aCoeff(0);
    if (std::abs(a0 - T(1)) < std::numeric_limits<T>::epsilon())
        return;

    m_aCoeff /= a0;
    m_bCoeff /= a0;
}

template <typename T>
void GenericFilter<T>::checkCoeffs(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff)
{
    Expects(aCoeff.size() > 0);
    Expects(std::abs(aCoeff[0]) > std::numeric_limits<T>::epsilon());
    Expects(bCoeff.size() > 0);
}

} // namespace difi