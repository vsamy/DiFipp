#define BOOST_TEST_MODULE GenericFilterTests

#include "fratio"
#include <boost/test/unit_test.hpp>
#include <vector>

BOOST_AUTO_TEST_CASE(FILTER_FAILURES)
{
    auto dfd = fratio::DigitalFilterd(Eigen::VectorXd(), Eigen::VectorXd::Constant(2, 0));
    BOOST_REQUIRE(dfd.status() == fratio::FilterStatus::A_COEFF_MISSING);
    dfd = fratio::DigitalFilterd(Eigen::VectorXd::Constant(2, 1), Eigen::VectorXd());
    BOOST_REQUIRE(dfd.status() == fratio::FilterStatus::B_COEFF_MISSING);
    dfd = fratio::DigitalFilterd(Eigen::VectorXd(), Eigen::VectorXd());
    BOOST_REQUIRE(dfd.status() == fratio::FilterStatus::ALL_COEFF_MISSING);
    dfd = fratio::DigitalFilterd(Eigen::VectorXd::Constant(2, 0), Eigen::VectorXd::Constant(2, 0));
    BOOST_REQUIRE(dfd.status() == fratio::FilterStatus::BAD_A_COEFF);
    dfd = fratio::DigitalFilterd();
    BOOST_REQUIRE(dfd.status() == fratio::FilterStatus::NONE);
    dfd = fratio::DigitalFilterd(Eigen::VectorXd::Constant(2, 1), Eigen::VectorXd::Constant(2, 0));
    BOOST_REQUIRE(dfd.status() == fratio::FilterStatus::READY);
    auto mad = fratio::MovingAveraged(0);
    BOOST_REQUIRE(mad.status() == fratio::FilterStatus::BAD_ORDER_SIZE);
    auto bfd = fratio::Butterworthd(0, 10, 100);
    BOOST_REQUIRE(bfd.status() == fratio::FilterStatus::BAD_ORDER_SIZE);
}
