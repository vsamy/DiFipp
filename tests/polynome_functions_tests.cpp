#define BOOST_TEST_MODULE polynome_functions_tests

#include "fratio.h"
#include "warning_macro.h"
#include <boost/test/unit_test.hpp>

struct SystemInt {
    std::vector<int> data = { 1, 1, 1, 1 };
    std::vector<int> results = { 1, -4, 6, -4, 1 };
};

struct SystemCInt {
    std::vector<std::complex<int>> data = { { 1, 1 }, { -1, 4 }, { 12, -3 }, { 5, 2 } };
    std::vector<std::complex<int>> results = { { 1, 0 }, { -17, -4 }, { 66, 97 }, { 127, -386 }, { -357, 153 } };
};

DISABLE_CONVERSION_WARNING_BEGIN

template <typename T>
struct SystemFloat {
    std::vector<T> data = { 0.32, -0.0518, 41.4, 0.89 };
    std::vector<T> results = { 1., -42.558199999999999, 48.171601999999993, -9.181098159999999, -0.610759296 };
};

template <typename T>
struct SystemCFloat {
    std::vector<std::complex<T>> data = { { 1.35, 0.2 }, { -1.5, 4.45 }, { 12.8, -3.36 }, { 5.156, 2.12 } };
    std::vector<std::complex<T>> results = { { 1., 0. }, { -17.806, -3.41 }, { 73.2776, 99.20074 }, { 101.857496, -444.634694 }, { -269.1458768, 388.7308864 } };
};

DISABLE_CONVERSION_WARNING_END

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_INT, SystemInt)
{
    auto res = fratio::VietaAlgoi::polyCoeffFromRoot(data);

    for (size_t i = 0; i < res.size(); ++i)
        BOOST_CHECK_EQUAL(res[i], results[i]);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_FLOAT, SystemFloat<float>)
{
    auto res = fratio::VietaAlgof::polyCoeffFromRoot(data);

    for (size_t i = 0; i < res.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(res[i] - results[i]), 1e-6f);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_DOUBLE, SystemFloat<double>)
{
    auto res = fratio::VietaAlgod::polyCoeffFromRoot(data);

    for (size_t i = 0; i < res.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(res[i] - results[i]), 1e-14);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CINT, SystemCInt)
{
    auto res = fratio::VietaAlgoci::polyCoeffFromRoot(data);

    for (size_t i = 0; i < res.size(); ++i)
        BOOST_CHECK_EQUAL(res[i], results[i]);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CFLOAT, SystemCFloat<float>)
{
    auto res = fratio::VietaAlgocf::polyCoeffFromRoot(data);

    for (size_t i = 0; i < res.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(res[i] - results[i]), 1e-4f);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CDOUBLE, SystemCFloat<double>)
{
    auto res = fratio::VietaAlgocd::polyCoeffFromRoot(data);

    for (size_t i = 0; i < res.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(res[i] - results[i]), 1e-12);
}
