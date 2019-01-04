// Copyright (c) 2019, Vincent SAMY
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the AIST nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#define BOOST_TEST_MODULE GenericFilterTests

#include "difi"
#include <boost/test/unit_test.hpp>
#include <vector>

BOOST_AUTO_TEST_CASE(FILTER_FAILURES)
{
    auto dfd = difi::DigitalFilterd(Eigen::VectorXd(), Eigen::VectorXd::Constant(2, 0));
    BOOST_REQUIRE(dfd.status() == difi::FilterStatus::A_COEFF_MISSING);
    dfd = difi::DigitalFilterd(Eigen::VectorXd::Constant(2, 1), Eigen::VectorXd());
    BOOST_REQUIRE(dfd.status() == difi::FilterStatus::B_COEFF_MISSING);
    dfd = difi::DigitalFilterd(Eigen::VectorXd(), Eigen::VectorXd());
    BOOST_REQUIRE(dfd.status() == difi::FilterStatus::ALL_COEFF_MISSING);
    dfd = difi::DigitalFilterd(Eigen::VectorXd::Constant(2, 0), Eigen::VectorXd::Constant(2, 0));
    BOOST_REQUIRE(dfd.status() == difi::FilterStatus::BAD_A_COEFF);
    dfd = difi::DigitalFilterd();
    BOOST_REQUIRE(dfd.status() == difi::FilterStatus::NONE);
    dfd = difi::DigitalFilterd(Eigen::VectorXd::Constant(2, 1), Eigen::VectorXd::Constant(2, 0));
    BOOST_REQUIRE(dfd.status() == difi::FilterStatus::READY);
    auto mad = difi::MovingAveraged(0);
    BOOST_REQUIRE(mad.status() == difi::FilterStatus::BAD_ORDER_SIZE);
    auto bfd = difi::Butterworthd(0, 10, 100);
    BOOST_REQUIRE(bfd.status() == difi::FilterStatus::BAD_ORDER_SIZE);
    bfd = difi::Butterworthd(5, 6, 5, 100);
    BOOST_REQUIRE(bfd.status() == difi::FilterStatus::BAD_BAND_FREQUENCY);
}
