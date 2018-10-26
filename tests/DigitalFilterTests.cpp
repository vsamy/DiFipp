#define BOOST_TEST_MODULE DigitalFilterTests

#include "fratio.h"
#include "warning_macro.h"
#include <boost/test/unit_test.hpp>

DISABLE_CONVERSION_WARNING_BEGIN

template <typename T>
struct System {
    std::vector<T> data = { 1, 2, 3, 4 };
    std::vector<T> aCoeff = { 1, -0.99993717 };
    std::vector<T> bCoeff = { 0.99996859, -0.99996859 };
    std::vector<T> results = { 0.99996859, 1.999874351973491, 2.999717289867956, 3.999497407630634 };
};

DISABLE_CONVERSION_WARNING_END

BOOST_FIXTURE_TEST_CASE(DIGITAL_FILTER_FLOAT, System<float>)
{
    auto df = fratio::DigitalFilterf(aCoeff, bCoeff);
    BOOST_REQUIRE_EQUAL(aCoeff.size(), df.aOrder());
    BOOST_REQUIRE_EQUAL(bCoeff.size(), df.bOrder());

    std::vector<float> filteredData;

    for (float d : data)
        filteredData.push_back(df.stepFilter(d));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-6f);

    df.resetFilter();
    filteredData = df.filter(data);
    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-6f);
}

BOOST_FIXTURE_TEST_CASE(DIGITAL_FILTER_DOUBLE, System<double>)
{
    auto df = fratio::DigitalFilterd(aCoeff, bCoeff);
    BOOST_REQUIRE_EQUAL(aCoeff.size(), df.aOrder());
    BOOST_REQUIRE_EQUAL(bCoeff.size(), df.bOrder());

    std::vector<double> filteredData;

    for (double d : data)
        filteredData.push_back(df.stepFilter(d));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-14);

    df.resetFilter();
    filteredData = df.filter(data);
    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-14);
}
