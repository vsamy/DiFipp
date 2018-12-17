#define BOOST_TEST_MODULE DigitalFilterTests

#include "fratio"
#include "test_functions.h"
#include "warning_macro.h"
#include <boost/test/unit_test.hpp>

DISABLE_CONVERSION_WARNING_BEGIN

template <typename T>
struct System {
    fratio::vectX_t<T> data = (fratio::vectX_t<T>(4) << 1, 2, 3, 4).finished();
    fratio::vectX_t<T> aCoeff = (fratio::vectX_t<T>(2) << 1, -0.99993717).finished();
    fratio::vectX_t<T> bCoeff = (fratio::vectX_t<T>(2) << 0.99996859, -0.99996859).finished();
    fratio::vectX_t<T> results = (fratio::vectX_t<T>(4) << 0.99996859, 1.999874351973491, 2.999717289867956, 3.999497407630634).finished();
};

DISABLE_CONVERSION_WARNING_END

BOOST_FIXTURE_TEST_CASE(DIGITAL_FILTER_FLOAT, System<float>)
{
    auto df = fratio::DigitalFilterf(aCoeff, bCoeff);
    test_coeffs(aCoeff, bCoeff, df);
    test_results(results, data, df);
}

BOOST_FIXTURE_TEST_CASE(DIGITAL_FILTER_DOUBLE, System<double>)
{
    auto df = fratio::DigitalFilterd(aCoeff, bCoeff);
    test_coeffs(aCoeff, bCoeff, df);
    test_results(results, data, df);
}
