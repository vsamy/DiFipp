#pragma once

#include "fratio"
#include <boost/test/unit_test.hpp>
#include <limits>

template <typename T>
void test_coeffs(const fratio::vectX_t<T>& aCoeff, const fratio::vectX_t<T>& bCoeff, const fratio::GenericFilter<T>& filter)
{
    BOOST_REQUIRE_EQUAL(aCoeff.size(), filter.aOrder());
    BOOST_REQUIRE_EQUAL(bCoeff.size(), filter.bOrder());
    fratio::vectX_t<T> faCoeff, fbCoeff;
    filter.getCoeffs(faCoeff, fbCoeff);
    for (Eigen::Index i = 0; i < faCoeff.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(aCoeff(i) - faCoeff(i)), std::numeric_limits<T>::epsilon() * 10);
    for (Eigen::Index i = 0; i < fbCoeff.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(bCoeff(i) - fbCoeff(i)), std::numeric_limits<T>::epsilon() * 10);
}

template <typename T>
void test_results(const fratio::vectX_t<T>& results, const fratio::vectX_t<T>& data, fratio::GenericFilter<T>& filter)
{
    fratio::vectX_t<T> filteredData(results.size());

    for (Eigen::Index i = 0; i < data.size(); ++i)
        filteredData(i) = filter.stepFilter(data(i));

    for (Eigen::Index i = 0; i < filteredData.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(filteredData(i) - results(i)), std::numeric_limits<T>::epsilon() * 1000);

    filter.resetFilter();
    filteredData = filter.filter(data);
    for (Eigen::Index i = 0; i < filteredData.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(filteredData(i) - results(i)), std::numeric_limits<T>::epsilon() * 1000);
}