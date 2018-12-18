#pragma once

#include "type_checks.h"
#include "typedefs.h"
#include <stddef.h>
#include <string>

namespace fratio {

template <typename T>
class GenericFilter {
    static_assert(std::is_floating_point<T>::value && !std::is_const<T>::value, "Only accept non-complex floating point types.");

public:
    static std::string filterStatus(FilterStatus status);

public:
    // Careful: Only an assert check for the filter status
    T stepFilter(const T& data);
    vectX_t<T> filter(const vectX_t<T>& data);
    bool getFilterResults(Eigen::Ref<vectX_t<T>> results, const vectX_t<T>& data);
    void resetFilter();

    template <typename T2>
    void setCoeffs(T2&& aCoeff, T2&& bCoeff);

    void getCoeffs(vectX_t<T>& aCoeff, vectX_t<T>& bCoeff) const;
    FilterStatus status() const noexcept { return m_status; }
    Eigen::Index aOrder() const noexcept { return m_aCoeff.size(); }
    Eigen::Index bOrder() const noexcept { return m_bCoeff.size(); }

protected:
    GenericFilter() = default;
    GenericFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff);
    virtual ~GenericFilter() = default;

    void normalizeCoeffs();
    bool checkCoeffs(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff);

protected:
    FilterStatus m_status;

private:
    vectX_t<T> m_aCoeff;
    vectX_t<T> m_bCoeff;
    vectX_t<T> m_filteredData;
    vectX_t<T> m_rawData;
};

} // namespace fratio

#include "GenericFilter.tpp"
