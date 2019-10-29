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

#include <limits>

namespace difi {

// Public functions

template <typename T>
void BaseFilter<T>::setType(Type type)
{
    Expects(type == Type::Centered ? m_bCoeff.size() > 2 && m_bCoeff.size() % 2 == 1 : true);
    m_center = (type == Type::Forward ? 0 : (m_bCoeff.size() - 1) / 2);
}

template <typename T>
template <typename T2>
void BaseFilter<T>::setCoeffs(T2&& aCoeff, T2&& bCoeff)
{
    static_assert(std::is_convertible_v<T2, vectX_t<T>>, "The coefficients types should be convertible to vectX_t<T>");

    Expects(checkCoeffs(aCoeff, bCoeff, (m_center == 0 ? Type::Forward : Type::Centered)));
    m_aCoeff = aCoeff;
    m_bCoeff = bCoeff;
    normalizeCoeffs();
    resetFilter();
    m_isInitialized = true;
}

template <typename T>
void BaseFilter<T>::getCoeffs(vectX_t<T>& aCoeff, vectX_t<T>& bCoeff) const noexcept
{
    aCoeff = m_aCoeff;
    bCoeff = m_bCoeff;
}

// Protected functions

template <typename T>
BaseFilter<T>::BaseFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff, FilterType type)
    : m_aCoeff(aCoeff)
    , m_bCoeff(bCoeff)
    , m_filteredData(aCoeff.size())
    , m_rawData(bCoeff.size())
{
    Expects(checkCoeffs(aCoeff, bCoeff, type));
    m_center = (type == Type::Forward ? 0 : (bCoeff.size() - 1) / 2);
    normalizeCoeffs();
    resetFilter();
    m_isInitialized = true;
}

template <typename T>
void BaseFilter<T>::normalizeCoeffs()
{
    T a0 = m_aCoeff(0);
    if (std::abs(a0 - T(1)) < std::numeric_limits<T>::epsilon())
        return;

    m_aCoeff /= a0;
    m_bCoeff /= a0;
}

template <typename T>
bool BaseFilter<T>::checkCoeffs(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff, FilterType type)
{
    bool centering = (type == Type::Centered ? (bCoeff.size() % 2 == 1) : true);
    return aCoeff.size() > 0 && std::abs(aCoeff[0]) > std::numeric_limits<T>::epsilon() && bCoeff.size() > 0 && centering;
}

} // namespace difi