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

#pragma once
#include "differentiators.h"
#include "noisy_function_generator.h"
#include <catch2/catch.hpp>
#include <limits>

constexpr const int STEPS = 200;
constexpr const int SIN_FREQUENCY = 60;

using namespace difi;

template <typename T>
struct TestFun {
    vectX_t<T> f;
    vectX_t<T> d;
    int center;
    double eps;
};

template <typename T, typename Filter>
struct TestRunner {
    void operator()(const TestFun<T>& tf, Filter& f)
    {
        for (int i = 0; i < 50; ++i)
            f.step(tf.f(i)); // First initialize, some steps

        for (int i = 20; i < STEPS; ++i)
            REQUIRE_SMALL(std::abs(f.step(tf.f(i)) - tf.d(i - tf.center)), tf.eps);
    }
};

template <typename T, size_t Ind, typename... Ts>
struct Tester {
    void operator()(const TestFun<T>& tf, std::tuple<Ts...> f)
    {
        TestRunner(tf, std::get<Ind>(f));
        Tester<T, Ind - 1, Ts...>(tf, f);
    }
};

template <typename T, typename... Ts>
struct Tester<T, 0, Ts...> {
    void operator()(const TestFun<T>& tf, std::tuple<Ts...> f) 
    {
        TestRunner(tf, std::get<0>(f));
    }
};

template <typename T, size_t N>
using centralList = std::tuple<CentralDiff<T, N>, LowNoiseLanczosDiff<T, N>, SuperLowNoiseLanczosDiff<T, N>, CenteredNoiseRobust2Diff<T, N>, CenteredNoiseRobust4Diff<T, N>>;

TEMPLATE_TEST_CASE("Sinus time-fixed central derivative", "[sin][central]", float, double)
{

    FunctionGenerator<TestType> fg = sinGenerator<TestType>(STEPS, SIN_FREQUENCY);

    auto list7 = centralList<TestType, 7>{};
    TestFun<TestType> list7Param = {std::get<0>(fg), std::get<2>(fg), 3, std::numeric_limits<TestType>::epsilon() * 100};
    TestFun<TestType> list7NoisyParam = {std::get<1>(fg), std::get<2>(fg), 3, std::numeric_limits<TestType>::epsilon() * 100};
    auto list9 = centralList<TestType, 9>{};
    TestFun<TestType> list9Param = {std::get<0>(fg), std::get<2>(fg), 4, std::numeric_limits<TestType>::epsilon() * 100};
    TestFun<TestType> list9NoisyParam = {std::get<1>(fg), std::get<2>(fg), 4, std::numeric_limits<TestType>::epsilon() * 100};

    // Check no noisy function
    

    // Check for noisy function

}
