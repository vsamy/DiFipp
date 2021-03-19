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

#include "diffTesters.h"
#include "doctest/doctest.h"
#include "doctest_helper.h"
#include "noisy_function_generator.h"
#include <limits>

constexpr const int STEPS = 500;
constexpr const int SIN_AMPLITUDE = 100;
constexpr const int SIN_FREQUENCY = 1;
template <typename T> constexpr const std::array<T, 5> POLY_4 = { 2, -3, 3, 2, -5 };
template <typename T> constexpr const std::array<T, 8> POLY_7 = { 10, -12, 5, 6, -9, -3, -2, 4 };

template <int N>
vectN_t<double, N> generateLNLCoeffs()
{
    static_assert(N > 2 && N % 2 == 1, "N must be odd");
    vectN_t<double, N> lnl;
    const int M = (N - 1) / 2;
    for (int i = 0; i < M + 1; ++i) {
        lnl(i) = static_cast<double>(M - i);
        lnl(M + i) = -static_cast<double>(i);
    }

    return lnl;
}

template <int N>
void checkCoeffs(const vectN_t<double, N>& coeffs, const vectN_t<double, N>& goodCoeffs)
{
    for (int i = 0; i < N; ++i)
        REQUIRE_SMALL(std::abs(coeffs(i) - goodCoeffs(i)), std::numeric_limits<double>::epsilon() * 5);
}

TEST_CASE("Coefficient calculation") // Some coeffs are computed, the rest are given
{
    // LNL coeffs
    checkCoeffs<5>(details::GetLNLCoeffs<double, 5>{}(), generateLNLCoeffs<5>() / 10.);
    checkCoeffs<7>(details::GetLNLCoeffs<double, 7>{}(), generateLNLCoeffs<7>() / 28.);
    checkCoeffs<9>(details::GetLNLCoeffs<double, 9>{}(), generateLNLCoeffs<9>() / 60.);
    checkCoeffs<11>(details::GetLNLCoeffs<double, 11>{}(), generateLNLCoeffs<11>() / 110.);

    // FNR coeffs
    checkCoeffs<3>(details::GetFNRCoeffs<double, 3>{}(), (vectN_t<double, 3>() << 1., 0., -1.).finished() / 2.);
    checkCoeffs<4>(details::GetFNRCoeffs<double, 4>{}(), (vectN_t<double, 4>() << 1., 1., -1., -1.).finished() / 4.);
    checkCoeffs<5>(details::GetFNRCoeffs<double, 5>{}(), (vectN_t<double, 5>() << 1., 2., 0., -2., -1.).finished() / 8.);
    checkCoeffs<6>(details::GetFNRCoeffs<double, 6>{}(), (vectN_t<double, 6>() << 1., 3., 2., -2., -3., -1.).finished() / 16.);
    checkCoeffs<7>(details::GetFNRCoeffs<double, 7>{}(), (vectN_t<double, 7>() << 1., 4., 5., 0., -5., -4., -1.).finished() / 32.);
    checkCoeffs<8>(details::GetFNRCoeffs<double, 8>{}(), (vectN_t<double, 8>() << 1., 5., 9., 5., -5., -9., -5., -1.).finished() / 64.);
    checkCoeffs<9>(details::GetFNRCoeffs<double, 9>{}(), (vectN_t<double, 9>() << 1., 6., 14., 14., 0., -14., -14., -6., -1.).finished() / 128.);
    checkCoeffs<10>(details::GetFNRCoeffs<double, 10>{}(), (vectN_t<double, 10>() << 1., 7., 20., 28., 14., -14., -28., -20., -7., -1.).finished() / 256.);
    checkCoeffs<11>(details::GetFNRCoeffs<double, 11>{}(), (vectN_t<double, 11>() << 1., 8., 27., 48., 42., 0., -42., -48., -27., -8., -1.).finished() / 512.);
}

TEST_CASE("Sinus time-fixed central derivative")
{
    double dt = 0.001;
    auto sg = sinGenerator<double>(STEPS, SIN_AMPLITUDE, SIN_FREQUENCY, dt);
    auto ct7 = tester::central_list<7>{};
    auto ct9 = tester::central_list<9>{};

    tester::set_time_steps(ct7, dt);
    tester::set_time_steps(ct9, dt);

    difi::vectX_t<double> eps{ 5 };
    eps << 1e-10, 1e-1, 1e-6, 1e-1, 1e-6; // Value checked with MATLAB

    tester::run_tests(ct7, std::get<0>(sg), std::get<1>(sg), eps);
    tester::run_tests(ct9, std::get<0>(sg), std::get<1>(sg), eps);
}

TEST_CASE("Polynome time-fixed central derivative")
{
    double dt = 0.001;
    {
        auto pg = polyGenerator<double>(STEPS, POLY_4<double>, dt);
        auto ct7 = tester::central_list<7>{};
        auto ct9 = tester::central_list<9>{};

        tester::set_time_steps(ct7, dt);
        tester::set_time_steps(ct9, dt);

        difi::vectX_t<double> eps{ 5 };
        eps << 1e-12, 1e-4, 1e-12, 1e-4, 1e-12; // Value checked with MATLAB

        tester::run_tests(ct7, std::get<0>(pg), std::get<1>(pg), eps);
        tester::run_tests(ct9, std::get<0>(pg), std::get<1>(pg), eps);
    }
    {
        auto pg = polyGenerator<double>(STEPS, POLY_7<double>, dt);
        auto ct7 = tester::central_list<7>{};
        auto ct9 = tester::central_list<9>{};

        tester::set_time_steps(ct7, dt);
        tester::set_time_steps(ct9, dt);

        difi::vectX_t<double> eps{ 5 };
        eps << 1e-11, 1e-3, 1e-9, 1e-4, 1e-9; // Value checked with MATLAB

        tester::run_tests(ct7, std::get<0>(pg), std::get<1>(pg), eps);
        tester::run_tests(ct9, std::get<0>(pg), std::get<1>(pg), eps);
    }
}

TEST_CASE("2nd order sinus time-fixed center derivative")
{
    double dt = 0.001;
    auto sg = sinGenerator<double>(STEPS, SIN_AMPLITUDE, SIN_FREQUENCY, dt);
    auto ct7 = std::tuple<CenteredDiffSecondOrderd<7>>{};
    auto ct9 = std::tuple<CenteredDiffSecondOrderd<9>>{};
    auto ct11 = std::tuple<CenteredDiffSecondOrderd<11>>{};

    tester::set_time_steps(ct7, dt);
    tester::set_time_steps(ct9, dt);
    tester::set_time_steps(ct11, dt);

    difi::vectX_t<double> eps{ 1 };
    eps << 2e-1;

    tester::run_tests(ct7, std::get<0>(sg), std::get<2>(sg), eps);
    tester::run_tests(ct9, std::get<0>(sg), std::get<2>(sg), eps);
    tester::run_tests(ct11, std::get<0>(sg), std::get<2>(sg), eps);
}

TEST_CASE("Sinus time-varying central derivative")
{
    double dt = 0.001;
    auto sg = tvSinGenerator<double>(STEPS, SIN_AMPLITUDE, SIN_FREQUENCY, dt);
    auto ct7 = tester::tv_central_list<7>{};
    auto ct9 = tester::tv_central_list<9>{};

    difi::vectX_t<double> eps{ 2 };
    eps << 1., 1.;

    tester::run_tests(ct7, std::get<0>(sg), std::get<1>(sg), std::get<2>(sg), eps);
    tester::run_tests(ct9, std::get<0>(sg), std::get<1>(sg), std::get<2>(sg), eps);
}

TEST_CASE("Polynome time-varying central derivative")
{
    double dt = 0.001;
    {
        auto pg = tvPolyGenerator<double>(STEPS, POLY_4<double>, dt);
        auto ct7 = tester::tv_central_list<7>{};
        auto ct9 = tester::tv_central_list<9>{};

        difi::vectX_t<double> eps{ 2 };
        eps << 1e-3, 1e-3;

        tester::run_tests(ct7, std::get<0>(pg), std::get<1>(pg), std::get<2>(pg), eps);
        tester::run_tests(ct9, std::get<0>(pg), std::get<1>(pg), std::get<2>(pg), eps);
    }
    {
        auto pg = tvPolyGenerator<double>(STEPS, POLY_7<double>, dt);
        auto ct7 = tester::tv_central_list<7>{};
        auto ct9 = tester::tv_central_list<9>{};

        difi::vectX_t<double> eps{ 2 };
        eps << 1e-2, 1e-2;

        tester::run_tests(ct7, std::get<0>(pg), std::get<1>(pg), std::get<2>(pg), eps);
        tester::run_tests(ct9, std::get<0>(pg), std::get<1>(pg), std::get<2>(pg), eps);
    }
}

// TEST_CASE("2nd order sinus time-varying center derivative", "[tv][sin][center][2nd]")
// {
//     // Test not passing.
//     double dt = 0.001;
//     auto sg = tvSinGenerator<double>(STEPS, SIN_AMPLITUDE, SIN_FREQUENCY, dt);
//     auto ct7 = generateTVC2OTester<double, 7>(std::get<0>(sg), std::get<1>(sg), std::get<3>(sg));
//     auto ct9 = generateTVC2OTester<double, 9>(std::get<0>(sg), std::get<1>(sg), std::get<3>(sg));
//     auto ct11 = generateTVC2OTester<double, 11>(std::get<0>(sg), std::get<1>(sg), std::get<3>(sg));

//     std::array<double, 1> eps = { 2e-1 };
//     ct7.run(eps);
//     ct9.run(eps);
//     ct11.run(eps);
// }
