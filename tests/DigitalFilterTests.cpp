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
#include "doctest/doctest.h"
#include "test_functions.h"
#include "warning_macro.h"

DISABLE_CONVERSION_WARNING_BEGIN

template <typename T>
struct System {
    difi::vectX_t<T> data = (difi::vectX_t<T>(4) << 1, 2, 3, 4).finished();
    difi::vectX_t<T> aCoeff = (difi::vectX_t<T>(2) << 1, -0.99993717).finished();
    difi::vectX_t<T> bCoeff = (difi::vectX_t<T>(2) << 0.99996859, -0.99996859).finished();
    difi::vectX_t<T> results = (difi::vectX_t<T>(4) << 0.99996859, 1.999874351973491, 2.999717289867956, 3.999497407630634).finished();
};

DISABLE_CONVERSION_WARNING_END

TEST_CASE_TEMPLATE("Digital filter", T, float, double)
{
    System<T> s;
    auto df = difi::DigitalFilter<T>(s.aCoeff, s.bCoeff);
    test_coeffs(s.aCoeff, s.bCoeff, df, std::numeric_limits<T>::epsilon() * 10);
    test_results(s.results, s.data, df, std::numeric_limits<T>::epsilon() * 10);
}
