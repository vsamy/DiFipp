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

#include "difi"
#include "doctest/doctest.h"
#include "doctest_helper.h"
#include <limits>

template <typename T>
void test_coeffs(const difi::vectX_t<T>& aCoeff, const difi::vectX_t<T>& bCoeff, const difi::GenericFilter<T>& filter, T prec)
{
    REQUIRE_EQUAL(aCoeff.size(), filter.aOrder());
    REQUIRE_EQUAL(bCoeff.size(), filter.bOrder());
    difi::vectX_t<T> faCoeff, fbCoeff;
    filter.getCoeffs(faCoeff, fbCoeff);
    for (Eigen::Index i = 0; i < faCoeff.size(); ++i)
        REQUIRE_SMALL(std::abs(aCoeff(i) - faCoeff(i)), prec);
    for (Eigen::Index i = 0; i < fbCoeff.size(); ++i)
        REQUIRE_SMALL(std::abs(bCoeff(i) - fbCoeff(i)), prec);
}

template <typename T>
void test_results(const difi::vectX_t<T>& results, const difi::vectX_t<T>& data, difi::GenericFilter<T>& filter, T prec)
{
    difi::vectX_t<T> filteredData(results.size());

    for (Eigen::Index i = 0; i < data.size(); ++i)
        filteredData(i) = filter.stepFilter(data(i));

    for (Eigen::Index i = 0; i < filteredData.size(); ++i)
        REQUIRE_SMALL(std::abs(filteredData(i) - results(i)), prec);

    filter.resetFilter();
    filteredData = filter.filter(data);
    for (Eigen::Index i = 0; i < filteredData.size(); ++i)
        REQUIRE_SMALL(std::abs(filteredData(i) - results(i)), prec);
}