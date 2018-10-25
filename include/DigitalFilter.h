#pragma once

#include "GenericFilter.h"

namespace fratio {

template <typename T>
class DigitalFilter : public GenericFilter<T> {
public:
    DigitalFilter() = default;
    DigitalFilter(const std::vector<T>& aCoeff, const std::vector<T>& bCoeff)
        : GenericFilter<T>(aCoeff, bCoeff)
    {
    }

    void setCoeff(const std::vector<double>& aCoeff, const std::vector<double>& bCoeff)
    {
        checkCoeff(aCoeff, bCoeff);

        m_aCoeff = aCoeff;
        m_bCoeff = bCoeff;
        resetFilter();

        normalize();
    }

    size_t aOrder() const noexcept { return m_aCoeff.size(); }
    size_t bOrder() const noexcept { return m_bCoeff.size(); }
};

} // namespace fratio