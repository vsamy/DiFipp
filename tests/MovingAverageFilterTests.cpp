#define BOOST_TEST_MODULE MovingAverageFilterTests

#include "fratio.h"
#include <boost/test/unit_test.hpp>

template <typename T>
struct System {
    Eigen::VectorX<T> data = (Eigen::VectorX<T>(6) << 1, 2, 3, 4, 5, 6).finished();
    size_t windowSize = 4;
    Eigen::VectorX<T> results = (Eigen::VectorX<T>(6) << 0.25, 0.75, 1.5, 2.5, 3.5, 4.5).finished();
};

BOOST_FIXTURE_TEST_CASE(MOVING_AVERAGE_FLOAT, System<float>)
{
    auto maf = fratio::MovingAveragef(windowSize);
    std::vector<float> filteredData;
    for (Eigen::Index i = 0; i < data.size(); ++i)
        filteredData.push_back(maf.stepFilter(data(i)));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-6f);

    maf.resetFilter();
    Eigen::VectorXf fData = maf.filter(data);
    for (Eigen::Index i = 0; i < fData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(fData(i) - results(i)), 1e-6f);
}

BOOST_FIXTURE_TEST_CASE(MOVING_AVERAGE_DOUBLE, System<double>)
{
    auto maf = fratio::MovingAveraged(windowSize);
    std::vector<double> filteredData;
    for (Eigen::Index i = 0; i < data.size(); ++i)
        filteredData.push_back(maf.stepFilter(data(i)));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-14);

    maf.resetFilter();
    Eigen::VectorXd fData = maf.filter(data);
    for (Eigen::Index i = 0; i < fData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(fData(i) - results(i)), 1e-14);
}
