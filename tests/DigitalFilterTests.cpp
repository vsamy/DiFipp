#define BOOST_TEST_MODULE DigitalFilterTests

#include "fratio.h"
#include "typedefs.h"
#include "warning_macro.h"
#include <boost/test/unit_test.hpp>

DISABLE_CONVERSION_WARNING_BEGIN

template <typename T>
struct System {
    Eigen::VectorX<T> data = (Eigen::VectorX<T>(4) << 1, 2, 3, 4).finished();
    Eigen::VectorX<T> aCoeff = (Eigen::VectorX<T>(4) << 1, -0.99993717).finished();
    Eigen::VectorX<T> bCoeff = (Eigen::VectorX<T>(4) << 0.99996859, -0.99996859).finished();
    Eigen::VectorX<T> results = (Eigen::VectorX<T>(4) << 0.99996859, 1.999874351973491, 2.999717289867956, 3.999497407630634).finished();
};

DISABLE_CONVERSION_WARNING_END

BOOST_FIXTURE_TEST_CASE(DIGITAL_FILTER_FLOAT, System<float>)
{
    auto df = fratio::DigitalFilterf(aCoeff, bCoeff);
    BOOST_REQUIRE_EQUAL(aCoeff.size(), df.aOrder());
    BOOST_REQUIRE_EQUAL(bCoeff.size(), df.bOrder());

    std::vector<float> filteredData;

    for (Eigen::Index i = 0; i < data.size(); ++i)
        filteredData.push_back(df.stepFilter(data(i)));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results(i)), 1e-6f);

    df.resetFilter();
    Eigen::VectorXf fData = df.filter(data);
    for (Eigen::Index i = 0; i < fData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(fData(i) - results(i)), 1e-6f);
}

BOOST_FIXTURE_TEST_CASE(DIGITAL_FILTER_DOUBLE, System<double>)
{
    auto df = fratio::DigitalFilterd(aCoeff, bCoeff);
    BOOST_REQUIRE_EQUAL(aCoeff.size(), df.aOrder());
    BOOST_REQUIRE_EQUAL(bCoeff.size(), df.bOrder());

    std::vector<double> filteredData;

    for (Eigen::Index i = 0; i < data.size(); ++i)
        filteredData.push_back(df.stepFilter(data(i)));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results(i)), 1e-14);

    df.resetFilter();
    Eigen::VectorXd fData = df.filter(data);
    for (Eigen::Index i = 0; i < fData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(fData(i) - results(i)), 1e-14);
}
