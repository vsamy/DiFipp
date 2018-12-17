#define BOOST_TEST_MODULE polynome_functions_tests

#include "fratio"
#include "warning_macro.h"
#include <boost/test/unit_test.hpp>

using c_int_t = std::complex<int>;
template <typename T>
using c_t = std::complex<T>;

struct SystemInt {
    Eigen::VectorX<int> data = (Eigen::VectorX<int>(4) << 1, 1, 1, 1).finished();
    Eigen::VectorX<int> results = (Eigen::VectorX<int>(5) << 1, -4, 6, -4, 1).finished();
};

struct SystemCInt {
    Eigen::VectorX<c_int_t> data = (Eigen::VectorX<c_int_t>(4) << c_int_t{ 1, 1 }, c_int_t{ -1, 4 }, c_int_t{ 12, -3 }, c_int_t{ 5, 2 }).finished();
    Eigen::VectorX<c_int_t> results = (Eigen::VectorX<c_int_t>(5) << c_int_t{ 1, 0 }, c_int_t{ -17, -4 }, c_int_t{ 66, 97 }, c_int_t{ 127, -386 }, c_int_t{ -357, 153 }).finished();
};

DISABLE_CONVERSION_WARNING_BEGIN

template <typename T>
struct SystemFloat {
    Eigen::VectorX<T> data = (Eigen::VectorX<T>(4) << 0.32, -0.0518, 41.4, 0.89).finished();
    Eigen::VectorX<T> results = (Eigen::VectorX<T>(5) << 1, -42.558199999999999, 48.171601999999993, -9.181098159999999, -0.610759296).finished();
};

template <typename T>
struct SystemCFloat {
    Eigen::VectorX<c_t<T>> data = (Eigen::VectorX<c_t<T>>(4) << c_t<T>{ 1.35, 0.2 }, c_t<T>{ -1.5, 4.45 }, c_t<T>{ 12.8, -3.36 }, c_t<T>{ 5.156, 2.12 }).finished();
    Eigen::VectorX<c_t<T>> results = (Eigen::VectorX<c_t<T>>(5) << c_t<T>{ 1, 0 }, c_t<T>{ -17.806, -3.41 }, c_t<T>{ 73.2776, 99.20074 }, c_t<T>{ 101.857496, -444.634694 }, c_t<T>{ -269.1458768, 388.7308864 }).finished();
};

DISABLE_CONVERSION_WARNING_END

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_INT, SystemInt)
{
    auto res = fratio::VietaAlgoi::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_CHECK_EQUAL(res(i), results(i));
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_FLOAT, SystemFloat<float>)
{
    auto res = fratio::VietaAlgof::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(res(i) - results(i)), 1e-6f);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_DOUBLE, SystemFloat<double>)
{
    auto res = fratio::VietaAlgod::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(res(i) - results(i)), 1e-14);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CINT, SystemCInt)
{
    auto res = fratio::VietaAlgoci::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_CHECK_EQUAL(res(i), results(i));
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CFLOAT, SystemCFloat<float>)
{
    auto res = fratio::VietaAlgocf::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(res(i) - results(i)), 1e-4f);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CDOUBLE, SystemCFloat<double>)
{
    auto res = fratio::VietaAlgocd::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(res(i) - results(i)), 1e-12);
}
