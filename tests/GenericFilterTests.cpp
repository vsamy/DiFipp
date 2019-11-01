// Copyright (c) 2019, Vincent SAMY
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies,
// either expressed or implied, of the FreeBSD Project.

#include "difi"
#include <exception>
#include <vector>
#include <catch2/catch.hpp>

TEST_CASE("Filter failures", "[fail]")
{
    // A coeff are missing
    REQUIRE_THROWS_AS(difi::DigitalFilterd(Eigen::VectorXd(), Eigen::VectorXd::Constant(2, 0)), std::logic_error);
    // B coeff are missing
    REQUIRE_THROWS_AS(difi::DigitalFilterd(Eigen::VectorXd::Constant(2, 1), Eigen::VectorXd()), std::logic_error);
    // aCoeff(0) = 0
    REQUIRE_THROWS_AS(difi::DigitalFilterd(Eigen::VectorXd::Constant(2, 0), Eigen::VectorXd::Constant(2, 0)), std::logic_error);
    // Filter left uninitialized
    REQUIRE_NOTHROW(difi::DigitalFilterd());
    auto df = difi::DigitalFilterd();
    // Filter data with uninitialized filter
    REQUIRE_THROWS_AS(df.stepFilter(10.), std::logic_error);
    // window <= 0
    REQUIRE_THROWS_AS(difi::MovingAveraged(0), std::logic_error);
    // order <= 0
    REQUIRE_THROWS_AS(difi::Butterworthd(0, 10, 100), std::logic_error);
    // fc > 2*fs
    REQUIRE_THROWS_AS(difi::Butterworthd(2, 60, 100), std::logic_error);
    // Upper frequency < lower frequency
    REQUIRE_THROWS_AS(difi::Butterworthd(2, 6, 5, 100), std::logic_error);

    // Ok
    REQUIRE_NOTHROW(difi::DigitalFilterd(Eigen::VectorXd::Constant(2, 1), Eigen::VectorXd::Constant(2, 0)));
}
