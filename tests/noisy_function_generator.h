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

#include "typedefs.h"
#include <array>
#include <cmath>
#include <tuple>

template <typename T>
using FunctionGenerator = std::tuple<difi::vectX_t<T>, difi::vectX_t<T>, difi::vectX_t<T>>;

template <typename T>
FunctionGenerator<T> sinGenerator(int nrSteps, T amplitude, T frequency, T dt)
{
    using namespace difi;

    vectX_t<T> f(nrSteps);
    vectX_t<T> df(nrSteps);
    vectX_t<T> ddf(nrSteps);

    for (int i = 0; i < nrSteps; ++i) {
        // function
        f(i) = amplitude * std::sin(2 * pi<T> * frequency * i * dt);

        // derivative 1
        df(i) = 2 * pi<T> * frequency * amplitude * std::cos(2 * pi<T> * frequency * i * dt);

        // derivative 2
        ddf(i) = -4 * pi<T> * pi<T> * frequency * frequency * amplitude * std::sin(2 * pi<T> * frequency * i * dt);
    }

    return { f, df, ddf };
}

template <typename T, size_t N>
FunctionGenerator<T> polyGenerator(int nrSteps, std::array<T, N> coeffs, T dt)
{
    using namespace difi;
    static_assert(N >= 2, "N must be superior to 20");

    auto computePoly = [](const auto& coeffs, T time) {
        auto recursiveComputation = [time, &coeffs](size_t i, T result, const auto& self) -> T {
            if (i > 0)
                return self(i - 1, time * result + coeffs[i - 1], self);
            else
                return result;
        };

        return recursiveComputation(coeffs.size(), 0, recursiveComputation);
    };

    std::array<T, N - 1> dCoeffs;
    for (size_t i = 1; i < N; ++i)
        dCoeffs[i - 1] = static_cast<T>(i) * coeffs[i];

    std::array<T, N - 2> ddCoeffs;
    for (size_t i = 1; i < N - 1; ++i)
        ddCoeffs[i - 1] = static_cast<T>(i) * dCoeffs[i];

    vectX_t<T> f(nrSteps);
    vectX_t<T> df(nrSteps);
    vectX_t<T> ddf(nrSteps);
    for (int i = 0; i < nrSteps; ++i) {
        // function
        f(i) = computePoly(coeffs, i * dt);

        // derivative 1
        df(i) = computePoly(dCoeffs, i * dt);

        // derivative 2
        ddf(i) = computePoly(ddCoeffs, i * dt);
    }

    return { f, df, ddf };
}

// TV sin generator

template <typename T>
using TVFunctionGenerator = std::tuple<difi::vectX_t<T>, difi::vectX_t<T>, difi::vectX_t<T>, difi::vectX_t<T>>;

template <typename T>
TVFunctionGenerator<T> tvSinGenerator(int nrSteps, T amplitude, T frequency, T meanDt)
{
    using namespace difi;

    vectX_t<T> t(nrSteps);
    vectX_t<T> f(nrSteps);
    vectX_t<T> df(nrSteps);
    vectX_t<T> ddf(nrSteps);

    t(0) = T(0);
    f(0) = T(0);
    df(0) = 2 * pi<T> * frequency * amplitude;
    ddf(0) = T(0);
    for (int i = 1; i < nrSteps; ++i) {
        // time
        if (i % 3 == 0)
            t(i) = t(i - 1) + 0.9 * meanDt;
        else if (i % 3 == 1)
            t(i) = t(i - 1) + meanDt;
        else
            t(i) = t(i - 1) + 1.1 * meanDt;

        // function
        f(i) = amplitude * std::sin(2 * pi<T> * frequency * t(i));

        // derivative 1
        df(i) = 2 * pi<T> * frequency * amplitude * std::cos(2 * pi<T> * frequency * t(i));

        // derivative 2
        ddf(i) = -4 * pi<T> * pi<T> * frequency * frequency * amplitude * std::sin(2 * pi<T> * frequency * t(i));
    }

    return { t, f, df, ddf };
}

template <typename T, size_t N>
TVFunctionGenerator<T> tvPolyGenerator(int nrSteps, std::array<T, N> coeffs, T meanDt)
{
    using namespace difi;
    static_assert(N >= 2, "N must be superior to 20");

    auto computePoly = [](const auto& coeffs, T time) {
        auto recursiveComputation = [time, &coeffs](size_t i, T result, const auto& self) -> T {
            if (i > 0)
                return self(i - 1, time * result + coeffs[i - 1], self);
            else
                return result;
        };

        return recursiveComputation(coeffs.size(), 0, recursiveComputation);
    };

    std::array<T, N - 1> dCoeffs;
    for (size_t i = 1; i < N; ++i)
        dCoeffs[i - 1] = static_cast<T>(i) * coeffs[i];

    std::array<T, N - 2> ddCoeffs;
    for (size_t i = 1; i < N - 1; ++i)
        ddCoeffs[i - 1] = static_cast<T>(i) * dCoeffs[i];

    vectX_t<T> t(nrSteps);
    vectX_t<T> f(nrSteps);
    vectX_t<T> df(nrSteps);
    vectX_t<T> ddf(nrSteps);
    t(0) = T(0);
    f(0) = computePoly(coeffs, t(0));
    df(0) = computePoly(dCoeffs, t(0));
    ddf(0) = computePoly(ddCoeffs, t(0));
    for (int i = 1; i < nrSteps; ++i) {
        // time
        if (i % 3 == 0)
            t(i) = t(i - 1) + 0.9 * meanDt;
        else if (i % 3 == 1)
            t(i) = t(i - 1) + meanDt;
        else
            t(i) = t(i - 1) + 1.1 * meanDt;

        // function
        f(i) = computePoly(coeffs, t(i));

        // derivative 1
        df(i) = computePoly(dCoeffs, t(i));

        // derivative 2
        ddf(i) = computePoly(ddCoeffs, t(i));
    }

    return { t, f, df, ddf };
}
