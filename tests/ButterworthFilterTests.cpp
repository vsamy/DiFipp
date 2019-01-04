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

#define BOOST_TEST_MODULE ButterworthFilterTests

// Note: In term of precision, LP > HP > BP ~= BR

#include "difi"
#include "test_functions.h"
#include "warning_macro.h"
#include <boost/test/unit_test.hpp>

DISABLE_CONVERSION_WARNING_BEGIN

template <typename T>
struct System {
    difi::vectX_t<T> data = (difi::vectX_t<T>(8) << 1, 2, 3, 4, 5, 6, 7, 8).finished();
    int order = 5;
    T fc = 10;
    T fs = 100;
    T fLower = 5;
    T fUpper = 15;
    // LP
    difi::vectX_t<T> lpACoeffRes = (difi::vectX_t<T>(6) << 1, -2.975422109745684, 3.806018119320413, -2.545252868330468, 0.881130075437837, -0.125430622155356).finished();
    difi::vectX_t<T> lpBCoeffRes = (difi::vectX_t<T>(6) << 0.001282581078961, 0.006412905394803, 0.012825810789607, 0.012825810789607, 0.006412905394803, 0.001282581078961).finished();
    difi::vectX_t<T> lpResults = (difi::vectX_t<T>(8) << 0.001282581078961, 0.012794287652606, 0.062686244350084, 0.203933712825708, 0.502244959135609, 1.010304217144175, 1.744652693589064, 2.678087381460197).finished();
    // HP
    difi::vectX_t<T> hpACoeffRes = (difi::vectX_t<T>(6) << 1, -2.975422109745683, 3.806018119320411, -2.545252868330467, 0.8811300754378368, -0.1254306221553557).finished();
    difi::vectX_t<T> hpBCoeffRes = (difi::vectX_t<T>(6) << 0.3541641810934298, -1.770820905467149, 3.541641810934299, -3.541641810934299, 1.770820905467149, -0.3541641810934298).finished();
    difi::vectX_t<T> hpResults = (difi::vectX_t<T>(8) << 0.3541641810934298, -0.008704608374924483, -0.3113626313910076, -0.3460321436983160, -0.1787600153274098, 0.04471440201428267, 0.2059279258827846, 0.2533941579793959).finished();
    // BP
    difi::vectX_t<T> bpACoeffRes = (difi::vectX_t<T>(11) << 1, -6.784299264603903, 21.577693329895588, -42.338550072279737, 56.729081385507655, -54.208087151300411, 37.399203252161037, -18.397491390111661, 6.180883710485754, -1.283022311577260, 0.125430622155356).finished();
    difi::vectX_t<T> bpBCoeffRes = (difi::vectX_t<T>(11) << 0.001282581078963, 0, -0.006412905394817, 0, 0.012825810789633, 0, -0.012825810789633, 0, 0.006412905394817, 0, -0.001282581078963).finished();
    difi::vectX_t<T> bpResults = (difi::vectX_t<T>(8) << 0.001282581078963, 0.011266576028733, 0.046195520115810, 0.116904647483408, 0.200574194600111, 0.232153315136604, 0.141350142008155, -0.086403129422609).finished();
    // BR
    difi::vectX_t<T> brACoeffRes = (difi::vectX_t<T>(11) << 1.000000000000000, -6.784299264603897, 21.577693329895553, -42.338550072279631, 56.729081385507484, -54.208087151300205, 37.399203252160873, -18.397491390111572, 6.180883710485723, -1.283022311577253, 0.125430622155356).finished();
    difi::vectX_t<T> brBCoeffRes = (difi::vectX_t<T>(11) << 0.354164181088899, -3.012700469326103, 12.021845263663774, -29.490886190815772, 49.130136704563000, -58.004276868015168, 49.130136704563000, -29.490886190815772, 12.021845263663774, -3.012700469326103, 0.354164181088899).finished();
    difi::vectX_t<T> brResults = (difi::vectX_t<T>(8) << 0.354164181088899, 0.098383686162154, 0.084355149987331, 0.375555141082278, 0.735622022349639, 1.008089442365644, 1.229578363722674, 1.537000959441760).finished();
};

DISABLE_CONVERSION_WARNING_END

BOOST_AUTO_TEST_CASE(FIND_BUTTERWORTH_LP_HP_FLOAT)
{
    // LP
    auto butterRequirement = difi::Butterworthf::findMinimumButter(40.f / 500.f, 150.f / 500.f, 3.f, 60.f);
    BOOST_REQUIRE_EQUAL(5, butterRequirement.first);
    BOOST_REQUIRE_SMALL(std::abs(static_cast<float>(0.081038494957764) - butterRequirement.second), std::numeric_limits<float>::epsilon() * 10);

    // HP
    butterRequirement = difi::Butterworthf::findMinimumButter(150.f / 500.f, 40.f / 500.f, 3.f, 60.f);
    BOOST_REQUIRE_EQUAL(5, butterRequirement.first);
    BOOST_REQUIRE_SMALL(std::abs(static_cast<float>(0.296655824107340) - butterRequirement.second), std::numeric_limits<float>::epsilon() * 10);
}

BOOST_AUTO_TEST_CASE(FIND_BUTTERWORTH_LP_HP_DOUBLE)
{
    // LP
    auto butterRequirement = difi::Butterworthd::findMinimumButter(40. / 500., 150. / 500., 3., 60.);
    BOOST_REQUIRE_EQUAL(5, butterRequirement.first);
    BOOST_REQUIRE_SMALL(std::abs(0.081038494957764 - butterRequirement.second), std::numeric_limits<double>::epsilon() * 10);

    // HP
    butterRequirement = difi::Butterworthd::findMinimumButter(150. / 500., 40. / 500., 3., 60.);
    BOOST_REQUIRE_EQUAL(5, butterRequirement.first);
    BOOST_REQUIRE_SMALL(std::abs(0.296655824107340 - butterRequirement.second), std::numeric_limits<double>::epsilon() * 10);
}

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_LP_FILTER_FLOAT, System<float>)
{
    auto bf = difi::Butterworthf(order, fc, fs);
    BOOST_REQUIRE_EQUAL(bf.aOrder(), bf.bOrder());
    test_coeffs(lpACoeffRes, lpBCoeffRes, bf, std::numeric_limits<float>::epsilon() * 10);
    test_results(lpResults, data, bf, std::numeric_limits<float>::epsilon() * 100);
}

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_LP_FILTER_DOUBLE, System<double>)
{
    auto bf = difi::Butterworthd(order, fc, fs);
    BOOST_REQUIRE_EQUAL(bf.aOrder(), bf.bOrder());
    test_coeffs(lpACoeffRes, lpBCoeffRes, bf, std::numeric_limits<double>::epsilon() * 10);
    test_results(lpResults, data, bf, std::numeric_limits<double>::epsilon() * 100);
}

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_HP_FILTER_FLOAT, System<float>)
{
    auto bf = difi::Butterworthf(order, fc, fs, difi::Butterworthf::Type::HighPass);
    BOOST_REQUIRE_EQUAL(bf.aOrder(), bf.bOrder());
    test_coeffs(hpACoeffRes, hpBCoeffRes, bf, std::numeric_limits<float>::epsilon() * 10);
    test_results(hpResults, data, bf, std::numeric_limits<float>::epsilon() * 1000);
}

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_HP_FILTER_DOUBLE, System<double>)
{
    auto bf = difi::Butterworthd(order, fc, fs, difi::Butterworthd::Type::HighPass);
    BOOST_REQUIRE_EQUAL(bf.aOrder(), bf.bOrder());
    test_coeffs(hpACoeffRes, hpBCoeffRes, bf, std::numeric_limits<double>::epsilon() * 10);
    test_results(hpResults, data, bf, std::numeric_limits<double>::epsilon() * 100);
}

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_BP_FILTER_FLOAT, System<float>)
{
    auto bf = difi::Butterworthf(order, fLower, fUpper, fs);
    BOOST_REQUIRE_EQUAL(bf.aOrder(), bf.bOrder());
    test_coeffs(bpACoeffRes, bpBCoeffRes, bf, std::numeric_limits<float>::epsilon() * 1000);
    test_results(bpResults, data, bf, std::numeric_limits<float>::epsilon() * 10000);
}

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_BP_FILTER_DOUBLE, System<double>)
{
    auto bf = difi::Butterworthd(order, fLower, fUpper, fs);
    BOOST_REQUIRE_EQUAL(bf.aOrder(), bf.bOrder());
    test_coeffs(bpACoeffRes, bpBCoeffRes, bf, std::numeric_limits<double>::epsilon() * 1000);
    test_results(bpResults, data, bf, std::numeric_limits<double>::epsilon() * 10000);
}

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_BR_FILTER_FLOAT, System<float>)
{
    auto bf = difi::Butterworthf(order, fLower, fUpper, fs, difi::Butterworthf::Type::BandReject);
    BOOST_REQUIRE_EQUAL(bf.aOrder(), bf.bOrder());
    test_coeffs(brACoeffRes, brBCoeffRes, bf, 1.f);
    test_results(brResults, data, bf, 1.f);
}

BOOST_FIXTURE_TEST_CASE(BUTTERWORTH_BR_FILTER_DOUBLE, System<double>)
{
    auto bf = difi::Butterworthd(order, fLower, fUpper, fs, difi::Butterworthd::Type::BandReject);
    BOOST_REQUIRE_EQUAL(bf.aOrder(), bf.bOrder());
    test_coeffs(brACoeffRes, brBCoeffRes, bf, 1e-8);
    test_results(brResults, data, bf, 1e-8);
}
