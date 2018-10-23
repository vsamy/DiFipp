// This file is part of copra.

// copra is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// copra is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with copra.  If not, see
// <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE FilterTest

#include "GenericFilter.h"
#include <Eigen/Core>
#include <boost/test/unit_test.hpp>

struct System1 {
    std::vector<double> data = { 1., 2., 3., 4. };
    std::vector<double> aCoeff = { 1., -0.99993717 };
    std::vector<double> bCoeff = { 0.99996859, -0.99996859 };
    std::vector<double> results = { 0.99996859, 1.999874351973491, 2.999717289867956, 3.999497407630634 };
};

BOOST_AUTO_TEST_CASE(FilterThrows)
{
    BOOST_REQUIRE_THROW(fratio::GenericFilter(1, {}, { 1., 2. }), std::runtime_error);
    BOOST_REQUIRE_THROW(fratio::GenericFilter(1, { 1., 2. }, {}), std::runtime_error);
    BOOST_REQUIRE_THROW(fratio::GenericFilter(1, { 0. }, { 1. }), std::invalid_argument);

    auto gf = fratio::GenericFilter();
    BOOST_REQUIRE_THROW(gf.setCoeff({}, { 1., 2. }), std::runtime_error);
    BOOST_REQUIRE_THROW(gf.setCoeff({ 1., 2. }, {}), std::runtime_error);
    BOOST_REQUIRE_THROW(gf.setCoeff({ 0. }, { 1. }), std::invalid_argument);
}

BOOST_FIXTURE_TEST_CASE(SimpleSystem, System1)
{
    auto gf = fratio::GenericFilter(1, aCoeff, bCoeff);
    BOOST_REQUIRE_EQUAL(aCoeff.size(), gf.aOrder());
    BOOST_REQUIRE_EQUAL(bCoeff.size(), gf.bOrder());

    std::vector<double> filteredData;

    for (double d : data)
        filteredData.push_back(gf.stepFilter(d));

    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-14);

    gf.resetFilter();
    filteredData = gf.filter(data);
    for (size_t i = 0; i < filteredData.size(); ++i)
        BOOST_CHECK_SMALL(std::abs(filteredData[i] - results[i]), 1e-14);
}
