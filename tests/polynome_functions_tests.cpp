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

#define BOOST_TEST_MODULE polynome_functions_tests

#include "difi"
#include "warning_macro.h"
#include <boost/test/unit_test.hpp>
#include <limits>

using c_int_t = std::complex<int>;
template <typename T>
using c_t = std::complex<T>;

struct SystemInt {
    difi::vectX_t<int> data = (difi::vectX_t<int>(4) << 1, 1, 1, 1).finished();
    difi::vectX_t<int> results = (difi::vectX_t<int>(5) << 1, -4, 6, -4, 1).finished();
};

struct SystemCInt {
    difi::vectX_t<c_int_t> data = (difi::vectX_t<c_int_t>(4) << c_int_t{ 1, 1 }, c_int_t{ -1, 4 }, c_int_t{ 12, -3 }, c_int_t{ 5, 2 }).finished();
    difi::vectX_t<c_int_t> results = (difi::vectX_t<c_int_t>(5) << c_int_t{ 1, 0 }, c_int_t{ -17, -4 }, c_int_t{ 66, 97 }, c_int_t{ 127, -386 }, c_int_t{ -357, 153 }).finished();
};

DISABLE_CONVERSION_WARNING_BEGIN

template <typename T>
struct SystemFloat {
    difi::vectX_t<T> data = (difi::vectX_t<T>(4) << 0.32, -0.0518, 41.4, 0.89).finished();
    difi::vectX_t<T> results = (difi::vectX_t<T>(5) << 1, -42.558199999999999, 48.171601999999993, -9.181098159999999, -0.610759296).finished();
};

template <typename T>
struct SystemCFloat {
    difi::vectX_t<c_t<T>> data = (difi::vectX_t<c_t<T>>(4) << c_t<T>{ 1.35, 0.2 }, c_t<T>{ -1.5, 4.45 }, c_t<T>{ 12.8, -3.36 }, c_t<T>{ 5.156, 2.12 }).finished();
    difi::vectX_t<c_t<T>> results = (difi::vectX_t<c_t<T>>(5) << c_t<T>{ 1, 0 }, c_t<T>{ -17.806, -3.41 }, c_t<T>{ 73.2776, 99.20074 }, c_t<T>{ 101.857496, -444.634694 }, c_t<T>{ -269.1458768, 388.7308864 }).finished();
};

DISABLE_CONVERSION_WARNING_END

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_INT, SystemInt)
{
    auto res = difi::VietaAlgoi::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_EQUAL(res(i), results(i));
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_FLOAT, SystemFloat<float>)
{
    auto res = difi::VietaAlgof::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(res(i) - results(i)), std::numeric_limits<float>::epsilon() * 1000);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_DOUBLE, SystemFloat<double>)
{
    auto res = difi::VietaAlgod::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(res(i) - results(i)), std::numeric_limits<double>::epsilon() * 1000);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CINT, SystemCInt)
{
    auto res = difi::VietaAlgoci::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_EQUAL(res(i), results(i));
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CFLOAT, SystemCFloat<float>)
{
    auto res = difi::VietaAlgocf::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(res(i) - results(i)), std::numeric_limits<float>::epsilon() * 1000);
}

BOOST_FIXTURE_TEST_CASE(POLYNOME_FUNCTION_CDOUBLE, SystemCFloat<double>)
{
    auto res = difi::VietaAlgocd::polyCoeffFromRoot(data);

    for (Eigen::Index i = 0; i < res.size(); ++i)
        BOOST_REQUIRE_SMALL(std::abs(res(i) - results(i)), std::numeric_limits<double>::epsilon() * 1000);
}
