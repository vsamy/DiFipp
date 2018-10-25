#pragma once

#include "type_checks.h"
#include <stddef.h>
#include <vector>

namespace fratio {

template <typename T, typename = std::enable_if_t<std::is_floating_point<T>::value && !std::is_const<T>::value>>
class GenericFilter;

template <typename T>
class GenericFilter<T> {
public:
    T stepFilter(T data);
    std::vector<T> filter(const std::vector<T>& data);
    void resetFilter();

    void getCoeff(std::vector<T>& aCoeff, std::vector<T>& bCoeff) const noexcept;

protected:
    GenericFilter() = default;
    GenericFilter(const std::vector<T>& aCoeff, const std::vector<T>& bCoeff);
    virtual ~GenericFilter() = default;

    void checkCoeff(const std::vector<T>& aCoeff, const std::vector<T>& bCoeff);
    void normalize();

protected:
    std::vector<T> m_aCoeff;
    std::vector<T> m_bCoeff;

private:
    std::vector<T> m_filteredData;
    std::vector<T> m_rawData;
};

} // namespace fratio

#include "GenericFilter.inl"