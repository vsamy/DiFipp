#define BOOST_TEST_MODULE MovingAverageFilterTests

#include "fratio"
#include "test_functions.h"
#include <boost/test/unit_test.hpp>

template <typename T>
struct System {
    fratio::vectX_t<T> data = (fratio::vectX_t<T>(6) << 1, 2, 3, 4, 5, 6).finished();
    int windowSize = 4;
    fratio::vectX_t<T> results = (fratio::vectX_t<T>(6) << 0.25, 0.75, 1.5, 2.5, 3.5, 4.5).finished();
};

BOOST_FIXTURE_TEST_CASE(MOVING_AVERAGE_FLOAT, System<float>)
{
    auto maf = fratio::MovingAveragef(windowSize);
    test_results(results, data, maf, std::numeric_limits<float>::epsilon() * 10);
}

BOOST_FIXTURE_TEST_CASE(MOVING_AVERAGE_DOUBLE, System<double>)
{
    auto maf = fratio::MovingAveraged(windowSize);
    test_results(results, data, maf, std::numeric_limits<double>::epsilon() * 10);
}
