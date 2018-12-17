#define BOOST_TEST_MODULE ButterworthFilterTests

#include "fratio"
#include "warning_macro.h"
#include <boost/test/unit_test.hpp>

DISABLE_CONVERSION_WARNING_BEGIN

template <typename T>
struct System {
    Eigen::VectorX<T> data = (Eigen::VectorX<T>(8) << 1, 2, 3, 4, 5, 6, 7, 8).finished();
    int order = 5;
    T fc = 10;
    T fs = 100;
    Eigen::VectorX<T> aCoeffRes = (Eigen::VectorX<T>(6) << 1.000000000000000, -2.975422109745684, 3.806018119320413, -2.545252868330468, 0.881130075437837, -0.125430622155356).finished();
    Eigen::VectorX<T> bCoeffRes = (Eigen::VectorX<T>(6) << 0.001282581078961, 0.006412905394803, 0.012825810789607, 0.012825810789607, 0.006412905394803, 0.001282581078961).finished();
    Eigen::VectorX<T> results = (Eigen::VectorX<T>(8) << 0.001282581078961, 0.012794287652606, 0.062686244350084, 0.203933712825708, 0.502244959135609, 1.010304217144175, 1.744652693589064, 2.678087381460197).finished();
};

DISABLE_CONVERSION_WARNING_END

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_FILTER_FLOAT, System<float>)
{
    auto bf = fratio::Butterworthf(order, fc, fs);

    std::vector<float> filteredData;
    Eigen::VectorX<float> aCoeff, bCoeff;
    bf.getCoeffs(aCoeff, bCoeff);

    BOOST_REQUIRE_EQUAL(aCoeff.size(), aCoeffRes.size());
    BOOST_REQUIRE_EQUAL(bCoeff.size(), bCoeffRes.size());
    BOOST_REQUIRE_EQUAL(aCoeff.size(), bCoeffRes.size());

    for (Eigen::Index i = 0; i < aCoeff.size(); ++i) {
        BOOST_CHECK_SMALL(std::abs(aCoeff(i) - aCoeffRes(i)), 1e-6f);
        BOOST_CHECK_SMALL(std::abs(bCoeff(i) - bCoeffRes(i)), 1e-6f);
    }

    for (Eigen::Index i = 0; i < data.size(); ++i)
        filteredData.push_back(bf.stepFilter(data(i)));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results(i)), 1e-6f);

    bf.resetFilter();
    Eigen::VectorXf fData = bf.filter(data);
    for (Eigen::Index i = 0; i < fData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(fData(i) - results(i)), 1e-6f);

    auto bf2 = fratio::Butterworthf();

    bf2.setFilterParameters(order, fc, fs);
    filteredData.clear();
    for (Eigen::Index i = 0; i < data.size(); ++i)
        filteredData.push_back(bf2.stepFilter(data(i)));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results(i)), 1e-6f);
}

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_FILTER_DOUBLE, System<double>)
{
    auto bf = fratio::Butterworthd(order, fc, fs);

    std::vector<double> filteredData;
    Eigen::VectorX<double> aCoeff, bCoeff;
    bf.getCoeffs(aCoeff, bCoeff);

    BOOST_REQUIRE_EQUAL(aCoeff.size(), aCoeffRes.size());
    BOOST_REQUIRE_EQUAL(bCoeff.size(), bCoeffRes.size());
    BOOST_REQUIRE_EQUAL(aCoeff.size(), bCoeffRes.size());

    for (Eigen::Index i = 0; i < aCoeff.size(); ++i) {
        BOOST_CHECK_SMALL(std::abs(aCoeff(i) - aCoeffRes(i)), 1e-14);
        BOOST_CHECK_SMALL(std::abs(bCoeff(i) - bCoeffRes(i)), 1e-14);
    }

    for (Eigen::Index i = 0; i < data.size(); ++i)
        filteredData.push_back(bf.stepFilter(data(i)));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results(i)), 1e-14);

    bf.resetFilter();
    Eigen::VectorXd fData = bf.filter(data);
    for (Eigen::Index i = 0; i < fData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(fData(i) - results(i)), 1e-14);

    auto bf2 = fratio::Butterworthd();

    bf2.setFilterParameters(order, fc, fs);
    filteredData.clear();
    for (Eigen::Index i = 0; i < data.size(); ++i)
        filteredData.push_back(bf2.stepFilter(data(i)));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results(i)), 1e-14);
}
