#define BOOST_TEST_MODULE GenericFilterTests

#include "fratio.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(FilterThrows)
{
    BOOST_REQUIRE_THROW(fratio::DigitalFilterd({}, { 1., 2. }), std::runtime_error);
    BOOST_REQUIRE_THROW(fratio::DigitalFilterd({ 1., 2. }, {}), std::runtime_error);
    BOOST_REQUIRE_THROW(fratio::DigitalFilterd({ 0. }, { 1. }), std::invalid_argument);

    auto df = fratio::DigitalFilterd();
    BOOST_REQUIRE_THROW(df.setCoeff({}, { 1., 2. }), std::runtime_error);
    BOOST_REQUIRE_THROW(df.setCoeff({ 1., 2. }, {}), std::runtime_error);
    BOOST_REQUIRE_THROW(df.setCoeff({ 0. }, { 1. }), std::invalid_argument);
}
