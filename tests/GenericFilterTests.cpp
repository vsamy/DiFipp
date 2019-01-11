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
#include <exception>
#include <vector>

BOOST_AUTO_TEST_CASE(FILTER_FAILURES)
{
    // A coeff are missing
    BOOST_REQUIRE_THROW(difi::DigitalFilterd(Eigen::VectorXd(), Eigen::VectorXd::Constant(2, 0)), std::logic_error);
    // B coeff are missing
    BOOST_REQUIRE_THROW(difi::DigitalFilterd(Eigen::VectorXd::Constant(2, 1), Eigen::VectorXd()), std::logic_error);
    // aCoeff(0) = 0
    BOOST_REQUIRE_THROW(difi::DigitalFilterd(Eigen::VectorXd::Constant(2, 0), Eigen::VectorXd::Constant(2, 0)), std::logic_error);
    // Filter left uninitialized
    BOOST_REQUIRE_NO_THROW(difi::DigitalFilterd());
    auto df = difi::DigitalFilterd();
    // Filter data with uninitialized filter
    BOOST_REQUIRE_THROW(df.stepFilter(10.), std::logic_error);
    // window <= 0
    BOOST_REQUIRE_THROW(difi::MovingAveraged(0), std::logic_error);
    // order <= 0
    BOOST_REQUIRE_THROW(difi::Butterworthd(0, 10, 100), std::logic_error);
    // fc > 2*fs
    BOOST_REQUIRE_THROW(difi::Butterworthd(2, 60, 100), std::logic_error);
    // Upper frequency < lower frequency
    BOOST_REQUIRE_THROW(difi::Butterworthd(2, 6, 5, 100), std::logic_error);

    // Ok
    BOOST_REQUIRE_NO_THROW(difi::DigitalFilterd(Eigen::VectorXd::Constant(2, 1), Eigen::VectorXd::Constant(2, 0)));
}
