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
#include "catch_helper.h"
#include "differentiators.h"
#include "difi"
#include <array>
#include <catch2/catch.hpp>

using namespace difi;

template <typename T, typename Tuple, typename Derived>
class DiffTester {
    friend Derived;
    static constexpr const size_t TSize = std::tuple_size_v<Tuple>;

public:
    DiffTester(const Tuple& filters, const vectX_t<T>& f, const vectX_t<T>& df)
        : m_filters(filters)
        , m_f(f)
        , m_df(df)
    {}

    void run(const std::array<T, TSize>& eps) { testRunner<TSize - 1>(eps); }

private:
    template <typename Filter>
    void test(Filter& filt, T eps)
    {
        static_cast<Derived&>(*this).testFilter(filt, eps);
    }

    template <size_t Index>
    void testRunner(const std::array<T, TSize>& eps)
    {
        test(std::get<Index>(m_filters), eps[Index]);
        testRunner<Index - 1>(eps);
    }

    template <>
    void testRunner<0>(const std::array<T, TSize>& eps) { test(std::get<0>(m_filters), eps[0]); }

private:
    vectX_t<T> m_f;
    vectX_t<T> m_df;
    Tuple m_filters;
};

template <typename T, typename Tuple>
class TFDiffTester : public DiffTester<T, Tuple, TFDiffTester<T, Tuple>> {
    friend DiffTester<T, Tuple, TFDiffTester<T, Tuple>>;

public:
    TFDiffTester(const Tuple& filters, const vectX_t<T>& f, const vectX_t<T>& df)
        : DiffTester(filters, f, df)
    {}

    void setFilterTimestep(T step) { setStep<TSize - 1>(step); }

private:
    template <typename Filter>
    void testFilter(Filter& filt, T eps)
    {
        for (int i = 0; i < 50; ++i) {
            auto value = filt.stepFilter(m_f(i)); // First initialize, some steps
            value = 2.;
        }

        for (int i = 50; i < STEPS; ++i) {
            auto value = filt.stepFilter(m_f(i));
            REQUIRE_SMALL(std::abs(value - m_df(i - filt.center())), eps);
        }
    }

    template <size_t Index>
    void setStep(T step)
    {
        std::get<Index>(m_filters).setTimestep(step);
        setStep<Index - 1>(step);
    }

    template <>
    void setStep<0>(T step) { std::get<0>(m_filters).setTimestep(step); }
};

template <typename T, typename Tuple>
class TVDiffTester : public DiffTester<T, Tuple, TVDiffTester<T, Tuple>> {
    friend DiffTester<T, Tuple, TVDiffTester<T, Tuple>>;

public:
    TVDiffTester(const Tuple& filters, const vectX_t<T>& t, const vectX_t<T>& f, const vectX_t<T>& df)
        : DiffTester(filters, f, df)
        , m_time(t)
    {}

private:
    template <typename Filter>
    void testFilter(Filter& filt, T eps)
    {
        for (int i = 0; i < 50; ++i) {
            auto value = filt.stepFilter(m_time(i), m_f(i)); // First initialize, some steps
            value = 2.;
        }

        for (int i = 50; i < STEPS; ++i) {
            auto value = filt.stepFilter(m_time(i), m_f(i));
            REQUIRE_SMALL(std::abs(value - m_df(i - filt.center())), eps);
        }
    }

private:
    vectX_t<T> m_time;
};

template <typename T, size_t N>
using central_list = std::tuple<CentralDiff<T, N>, LowNoiseLanczosDiff<T, N>, SuperLowNoiseLanczosDiff<T, N>, CenteredNoiseRobust2Diff<T, N>, CenteredNoiseRobust4Diff<T, N>>;

template <typename T, size_t N>
TFDiffTester<T, central_list<T, N>> generateCTester(const vectX_t<T>& f, const vectX_t<T>& df)
{
    return { central_list<T, N>{}, f, df };
}

template <typename T, size_t N>
TFDiffTester<T, std::tuple<CenteredSecondOrderDiff<T, N>>> generateC2OTester(const vectX_t<T>& f, const vectX_t<T>& df)
{
    return { std::tuple<CenteredSecondOrderDiff<T, N>>{}, f, df };
}

template <typename T, size_t N>
using tv_central_list = std::tuple<TVCenteredNoiseRobust2Diff<T, N>, TVCenteredNoiseRobust4Diff<T, N>>;

template <typename T, size_t N>
TVDiffTester<T, tv_central_list<T, N>> generateTVCTester(const vectX_t<T>& t, const vectX_t<T>& f, const vectX_t<T>& df)
{
    return { tv_central_list<T, N>{}, t, f, df };
}

template <typename T, size_t N>
TVDiffTester<T, std::tuple<TVCenteredSecondOrderDiff<T, N>>> generateTVC2OTester(const vectX_t<T>& t, const vectX_t<T>& f, const vectX_t<T>& df)
{
    return { std::tuple<TVCenteredSecondOrderDiff<T, N>>{}, t, f, df };
}