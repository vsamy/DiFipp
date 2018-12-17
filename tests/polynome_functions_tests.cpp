#define BOOST_TEST_MODULE polynome_functions_tests

#include "fratio"
#include "warning_macro.h"
#include <boost/test/unit_test.hpp>
#include <limits>

using c_int_t = std::complex<int>;
template <typename T>
using c_t = std::complex<T>;

struct SystemInt {
    fratio::vectX_t<int> data = (fratio::vectX_t<int>(4) << 1, 1, 1, 1).finished();
    fratio::vectX_t<int> results = (fratio::vectX_t<int>(5) << 1, -4, 6, -4, 1).finished();
};

struct SystemCInt {
    fratio::vectX_t<c_int_t> data = (fratio::vectX_t<c_int_t>(4) << c_int_t{ 1, 1 }, c_int_t{ -1, 4 }, c_int_t{ 12, -3 }, c_int_t{ 5, 2 }).finished();
    fratio::vectX_t<c_int_t> results = (fratio::vectX_t<c_int_t>(5) << c_int_t{ 1, 0 }, c_int_t{ -17, -4 }, c_int_t{ 66, 97 }, c_int_t{ 127, -386 }, c_int_t{ -357, 153 }).finished();
};

DISABLE_CONVERSION_WARNING_BEGIN

template <typename T>
struct SystemFloat {
    fratio::vectX_t<T> data = (fratio::vectX_t<T>(4) << 0.32, -0.0518, 41.4, 0.89).finished();
    fratio::vectX_t<T> results = (fratio::vectX_t<T>(5) << 1, -42.558199999999999, 48.171601999999993, -9.181098159999999, -0.610759296).finished();
};

template <typename T>
struct SystemCFloat {
    fratio::vectX_t<c_t<T>> data = (fratio::vectX_t<c_t<T>>(4) << c_t<T>{ 1.35, 0.2 }, c_t<T>{ -1.5, 4.45 }, c_t<T>{ 12.8, -3.36 }, c_t<T>{ 5.156, 2.12 }).finished();
    fratio::vectX_t<c_t<T>> results = (fratio::vectX_t<c_t<T>>(5) << c_t<T>{ 1, 0 }, c_t<T>{ -17.806, -3.41 }, c_t<T>{ 73.2776, 99.20074 }, c_t<T>{ 101.857496, -444.634694 }, c_t<T>{ -269.1458768, 388.7308864 }).finished();
};

DISABLE_CONVERSION_WARNING_END

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_INT, SystemInt)
{
    auto res = fratio::VietaAlgoi::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_EQUAL(res(i), results(i));
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_FLOAT, SystemFloat<float>)
{
    auto res = fratio::VietaAlgof::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(res(i) - results(i)), std::numeric_limits<float>::epsilon() * 1000);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_DOUBLE, SystemFloat<double>)
{
    auto res = fratio::VietaAlgod::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(res(i) - results(i)), std::numeric_limits<double>::epsilon() * 1000);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CINT, SystemCInt)
{
    auto res = fratio::VietaAlgoci::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_EQUAL(res(i), results(i));
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CFLOAT, SystemCFloat<float>)
{
    auto res = fratio::VietaAlgocf::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(res(i) - results(i)), std::numeric_limits<float>::epsilon() * 1000);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CDOUBLE, SystemCFloat<double>)
{
    auto res = fratio::VietaAlgocd::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(res(i) - results(i)), std::numeric_limits<double>::epsilon() * 1000);
}
