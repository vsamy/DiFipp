#pragma once

#include <complex>

namespace fratio {

template <typename T>
struct is_complex_t : public std::false_type {
};

template <typename T>
struct is_complex_t<std::complex<T>> : public std::true_type {
};

} // namespace fratio