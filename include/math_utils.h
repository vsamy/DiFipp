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
#include <type_traits>

namespace difi {

// https://stackoverflow.com/questions/44718971/calculate-binomial-coffeficient-very-reliably
template <typename T>
constexpr T Binomial(int n, int k)
{
    if (k > n)
        return T(0);
    else if (k == 0 || k == n)
        return T(1);
    else if (k == 1 || k == n - 1)
        return T(n);
    else if (2 * k < n)
        return Binomial<T>(n - 1, k - 1) * n / k;
    else
        return Binomial<T>(n - 1, k) * n / (n - k);
}

template <typename T>
constexpr T pow(T n, T k)
{
    static_assert(std::is_integral<T>::value, "For integral values only, eitherwise, use std::pow");
    if (n == 0)
        return (k == 0) ? 1 : 0;
    else if (k == 0)
        return T(1);
    else
        return n * pow(n, k - 1);
}

} // namespace difi
