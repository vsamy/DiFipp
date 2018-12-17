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
    Eigen::VectorX<T> filter(const Eigen::VectorX<T>& data);
    bool getFilterResults(Eigen::Ref<Eigen::VectorX<T>> results, const Eigen::VectorX<T>& data);
    void resetFilter();

    template <typename T2>
    bool setCoeffs(T2&& aCoeff, T2&& bCoeff);

    void getCoeffs(Eigen::VectorX<T>& aCoeff, Eigen::VectorX<T>& bCoeff) const;
    FilterStatus status() const noexcept { return m_status; }
    Eigen::Index aOrder() const noexcept { return m_aCoeff.size(); }
    Eigen::Index bOrder() const noexcept { return m_bCoeff.size(); }

protected:
    GenericFilter() = default;
    GenericFilter(const Eigen::VectorX<T>& aCoeff, const Eigen::VectorX<T>& bCoeff);
    virtual ~GenericFilter() = default;

    void normalizeCoeffs();
    bool checkCoeffs(const Eigen::VectorX<T>& aCoeff, const Eigen::VectorX<T>& bCoeff);

protected:
    FilterStatus m_status;

private:
    Eigen::VectorX<T> m_aCoeff;
    Eigen::VectorX<T> m_bCoeff;
    Eigen::VectorX<T> m_filteredData;
    Eigen::VectorX<T> m_rawData;
};

} // namespace fratio

#include "GenericFilter.tpp"
