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
#include <tuple>
#include <cmath>
#include <random>

template <typename T>
using FunctionGenerator = std::tuple<difi::vectX_t<T>, difi::vectX_t<T>, difi::vectX_t<T>>;

template <typename T>
FunctionGenerator<T> sinGenerator(int nrSteps, T frequency, T dt = 0.001)
{
    using namespace difi;

    std::random_device rd{};
    std::mt19937 gen{rd()};

    vectX_t<T> truth;
    vectX_t<T> noisy;
    vectX_t<T> derivative;

    for (int i = 0; i < nrSteps; ++i) {
        // truth
        truth(i) = std::sin(2 * pi<T> * frequency * i * dt);

        // noisy
        std::normal_distribution<T> d{truth(i), T(0.01)};
        noisy(i) = truth(i) + d(gen);

        // derivative
        derivative(i) = 2 * pi<T> * frequency * i * std::cos(2 * pi<T> * frequency * i * dt);
    }

    return { truth, noisy, derivative };
}

template <typename T>
FunctionGenerator<T> polyGenerator(int nrSteps, difi::vectX_t<T> coeffs, T dt = 0.001)
{
    using namespace difi;
    Expects(coeffs.size() >=2);
    std::random_device rd{};
    std::mt19937 gen{rd()};

    auto computePoly = [](const VectX_t<T>& coeffs, T time) {
        auto recursiveComputation = [time, &coeffs](int i, T result) {
            if (i > 0)
                return recursiveComputation(i - 1, time * result + coeffs(i - 1));
            else 
                return result;
        };

        return recursiveComputation(coeffs.size(), 0);
    };

    vectX_t<T> derivativeCoeffs(coeffs.size() - 1);
    for (Eigen::Index i = 1; i < coeffs.size(); ++i)
        derivativeCoeffs(i - 1) = i * coeffs.tail(i);

    vectX_t<T> truth;
    vectX_t<T> noisy;
    vectX_t<T> derivative;
    for (int i = 0; i < nrSteps; ++i) {
        // truth
        truth(i) = computePoly(coeffs, i * dt);

        // noisy
        std::normal_distribution<T> d{truth(i), T(0.01)};
        noisy(i) = truth(i) + d(gen);

        // derivative
        derivative(i) = computePoly(derivativeCoeffs, i * dt);
    }

    return { truth, noisy, derivative };
}