#define BOOST_TEST_MODULE ButterworthFilterTests

#include "fratio.h"
#include "warning_macro.h"
#include <boost/test/unit_test.hpp>

DISABLE_CONVERSION_WARNING_BEGIN

template <typename T>
struct System {
    std::vector<T> data = { 1., 2., 3., 4., 5., 6., 7., 8. };
    size_t order = 5;
    T fc = 10;
    T fs = 100;
    std::vector<T> aCoeffRes = { 1.000000000000000, -2.975422109745684, 3.806018119320413, -2.545252868330468, 0.881130075437837, -0.125430622155356 };
    std::vector<T> bCoeffRes = { 0.001282581078961, 0.006412905394803, 0.012825810789607, 0.012825810789607, 0.006412905394803, 0.001282581078961 };
    std::vector<T> results = { 0.001282581078961, 0.012794287652606, 0.062686244350084, 0.203933712825708, 0.502244959135609, 1.010304217144175, 1.744652693589064, 2.678087381460197 };
};

DISABLE_CONVERSION_WARNING_END

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_FILTER_FLOAT, System<float>)
{
    auto gf = fratio::Butterworthf(order, fc, fs);

    std::vector<float> aCoeff, bCoeff, filteredData;
    gf.getCoeff(aCoeff, bCoeff);

    BOOST_REQUIRE_EQUAL(aCoeff.size(), aCoeffRes.size());
    BOOST_REQUIRE_EQUAL(bCoeff.size(), bCoeffRes.size());
    BOOST_REQUIRE_EQUAL(aCoeff.size(), bCoeffRes.size());

    for (size_t i = 0; i < aCoeff.size(); ++i) {
        BOOST_CHECK_SMALL(std::abs(aCoeff[i] - aCoeffRes[i]), 1e-4f);
        BOOST_CHECK_SMALL(std::abs(bCoeff[i] - bCoeffRes[i]), 1e-4f);
    }

    for (float d : data)
        filteredData.push_back(gf.stepFilter(d));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-4f);

    gf.resetFilter();
    filteredData = gf.filter(data);
    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-4f);
}

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_FILTER_DOUBLE, System<double>)
{
    auto gf = fratio::Butterworthd(order, fc, fs);

    std::vector<double> aCoeff, bCoeff, filteredData;
    gf.getCoeff(aCoeff, bCoeff);

    BOOST_REQUIRE_EQUAL(aCoeff.size(), aCoeffRes.size());
    BOOST_REQUIRE_EQUAL(bCoeff.size(), bCoeffRes.size());
    BOOST_REQUIRE_EQUAL(aCoeff.size(), bCoeffRes.size());

    for (size_t i = 0; i < aCoeff.size(); ++i) {
        BOOST_CHECK_SMALL(std::abs(aCoeff[i] - aCoeffRes[i]), 1e-14);
        BOOST_CHECK_SMALL(std::abs(bCoeff[i] - bCoeffRes[i]), 1e-14);
    }

    for (double d : data)
        filteredData.push_back(gf.stepFilter(d));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-14);

    gf.resetFilter();
    filteredData = gf.filter(data);
    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-14);
}
