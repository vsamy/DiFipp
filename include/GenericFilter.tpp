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

namespace difi {

template <typename T>
T GenericFilter<T>::stepFilter(const T& data)
{
    Expects(m_isInitialized);

    // Slide data (can't use SIMD, but should be small)
    for (Eigen::Index i = m_rawData.size() - 1; i > 0; --i)
        m_rawData(i) = m_rawData(i - 1);
    for (Eigen::Index i = m_filteredData.size() - 1; i > 0; --i)
        m_filteredData(i) = m_filteredData(i - 1);

    m_rawData(0) = data;
    m_filteredData(0) = 0;
    m_filteredData(0) = m_bCoeff.dot(m_rawData) - m_aCoeff.dot(m_filteredData);
    return m_filteredData(0);
}

template <typename T>
vectX_t<T> GenericFilter<T>::filter(const vectX_t<T>& data)
{
    Expects(m_isInitialized);
    vectX_t<T> results(data.size());
    for (Eigen::Index i = 0; i < data.size(); ++i)
        results(i) = stepFilter(data(i));
    return results;
}

template <typename T>
void GenericFilter<T>::resetFilter() noexcept
{
    m_filteredData.setZero(m_aCoeff.size());
    m_rawData.setZero(m_bCoeff.size());
}

template <typename T>
T TVGenericFilter<T>::stepFilter(const T& time, const T& data)
{
    Expects(m_isInitialized);

    // Slide data (can't use SIMD, but should be small)
    for (Eigen::Index i = m_rawData.size() - 1; i > 0; --i) {
        m_rawData(i) = m_rawData(i - 1);
        m_timers(i) = m_timers(i - 1);
    }
    for (Eigen::Index i = m_filteredData.size() - 1; i > 0; --i)
        m_filteredData(i) = m_filteredData(i - 1);

    m_timers(0) = time;
    const Eigen::Index M = (m_rawData.size() - 1) / 2;
    auto bCoeff = m_bCoeff;
    for (Eigen::Index i = 1; i < M + 1; ++i) {
        const T diff = std::pow(m_timers(M - i) - m_timers(M + i), m_diffOrder);
        bCoeff(M + i) /= diff;
        bCoeff(M - i) /= diff;
        bCoeff(M) -= (bCoeff(M - i) + bCoeff(M + i));
    }
    m_rawData(0) = data;
    m_filteredData(0) = 0;
    m_filteredData(0) = bCoeff.dot(m_rawData) - m_aCoeff.dot(m_filteredData);
    return m_filteredData(0);
}

template <typename T>
vectX_t<T> TVGenericFilter<T>::filter(const vectX_t<T>& data, const vectX_t<T>& time)
{
    Expects(m_isInitialized);
    Expects(data.size() == time.size());
    vectX_t<T> results(data.size());
    for (Eigen::Index i = 0; i < data.size(); ++i)
        results(i) = stepFilter(time(i), data(i));
    return results;
}

template <typename T>
void TVGenericFilter<T>::resetFilter() noexcept
{
    m_filteredData.setZero(m_aCoeff.size());
    m_rawData.setZero(m_bCoeff.size());
    m_timers.setZero(m_bCoeff.size());
}

} // namespace difi