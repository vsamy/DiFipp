#define BOOST_TEST_MODULE MovingAverageFilterTests

#include "fratio.h"
#include <boost/test/unit_test.hpp>

template <typename T>
struct System {
    std::vector<T> data = { 1., 2., 3., 4., 5., 6. };
    size_t windowSize = 4;
    std::vector<T> results = { 0.25, 0.75, 1.5, 2.5, 3.5, 4.5 };
};

BOOST_FIXTURE_TEST_CASE(MOVING_AVERAGE_FLOAT, System<float>)
{
    auto ma = fratio::MovingAveragef(windowSize);
    std::vector<float> filteredData;
    for (float d : data)
        filteredData.push_back(ma.stepFilter(d));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-6f);

    ma.resetFilter();
    filteredData = ma.filter(data);
    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-6f);
}

BOOST_FIXTURE_TEST_CASE(MOVING_AVERAGE_DOUBLE, System<double>)
{
    auto ma = fratio::MovingAveraged(windowSize);
    std::vector<double> filteredData;
    for (double d : data)
        filteredData.push_back(ma.stepFilter(d));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-14);

    ma.resetFilter();
    filteredData = ma.filter(data);
    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-14);
}
