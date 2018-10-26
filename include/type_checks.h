#pragma once

#include <complex>

namespace fratio {

namespace internal {

    // Based on https://stackoverflow.com/a/41449096/9774052
    template <typename T>
    std::true_type is_tt_impl(std::complex<T>);
    std::false_type is_tt_impl(...);

    template <typename T>
    using is_complex = decltype(is_tt_impl(std::declval<T>()));

    template <typename T, bool B>
    struct sub_type {
        using type = T;
    };

    template <typename T>
    struct sub_type<T, true> {
        using type = typename T::value_type;
    };

    template <typename T>
    using complex_sub_type_t = typename sub_type<T, is_complex<T>::value>::type;

} // namespace internal

} // namespace fratio