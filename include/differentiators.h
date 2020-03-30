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
#include "DigitalFilter.h"
#include "math_utils.h"

// Read this if you are an adulator of the math god: https://arxiv.org/pdf/1709.08321.pdf

namespace difi {

namespace details {

// T: type
// N: Number of points

// Centered differentiators: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
template <typename T, int N> struct GetCDCoeffs {
    vectN_t<T, N> operator()() const;
};
template <typename T> struct GetCDCoeffs<T, 3> {
    vectN_t<T, 3> operator()() const { return (vectN_t<T, 3>() << T(1), T(0), T(-1)).finished() / T(2); }
};
template <typename T> struct GetCDCoeffs<T, 5> {
    vectN_t<T, 5> operator()() const { return (vectN_t<T, 5>() << T(-1), T(8), T(0), T(8), T(1)).finished() / T(12); }
};
template <typename T> struct GetCDCoeffs<T, 7> {
    vectN_t<T, 7> operator()() const { return (vectN_t<T, 7>() << T(1), T(-9), T(45), T(0), T(-45), T(9), T(-1)).finished() / T(60); }
};
template <typename T> struct GetCDCoeffs<T, 9> {
    vectN_t<T, 9> operator()() const { return (vectN_t<T, 9>() << T(-3), T(32), T(-168), T(672), T(0), T(-672), T(168), T(-32), T(3)).finished() / T(840); }
};

// Low-noise Lanczos differentiators: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/lanczos-low-noise-differentiators/
template <typename T, int N>
struct GetLNLCoeffs {
    vectN_t<T, N> operator()() const
    {
        static_assert(N > 2 && N % 2 == 1, "'N' must be odd.");
        constexpr const int M = (N - 1) / 2;
        constexpr const int Den = M * (M + 1) * (2 * M + 1);

        vectN_t<T, N> v{};
        v(M) = T(0);
        for (int k = 0; k < M; ++k) {
            v(k) = T(3) * static_cast<T>(M - k) / static_cast<T>(Den);
            v(N - k - 1) = -v(k);
        }
        return v;
    }
};

// Super Low-noise Lanczos differentiators: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/lanczos-low-noise-differentiators/
template <typename T, int N> struct GetSLNLCoeffs {
    vectN_t<T, N> operator()() const;
};
template <typename T> struct GetSLNLCoeffs<T, 7> {
    vectN_t<T, 7> operator()() const { return (vectN_t<T, 7>() << T(-22), T(67), T(58), T(0), T(-58), T(-67), T(22)).finished() / T(252); }
};
template <typename T> struct GetSLNLCoeffs<T, 9> {
    vectN_t<T, 9> operator()() const { return (vectN_t<T, 9>() << T(-86), T(142), T(193), T(126), T(0), T(-126), T(-193), T(-142), T(86)).finished() / T(1188); }
};
template <typename T> struct GetSLNLCoeffs<T, 11> {
    vectN_t<T, 11> operator()() const { return (vectN_t<T, 11>() << T(-300), T(294), T(532), T(503), T(296), T(0), T(-296), T(-503), T(-532), T(-294), T(300)).finished() / T(5148); }
};

// Backward Noise-Robust differentiators; http://www.holoborodko.com/pavel/wp-content/uploads/OneSidedNoiseRobustDifferentiators.pdf
template <typename T, int N>
struct GetFNRCoeffs {
    vectN_t<T, N> operator()() const
    {
        static_assert(N >= 2, "N should be greater than 2");
        constexpr const int BinCoeff = N - 2;
        constexpr const int M = N / 2;

        vectN_t<T, N> v{};
        v(0) = T(1);
        v(M) = T(0);
        v(N - 1) = T(-1);
        for (int i = 1; i < M; ++i) {
            v(i) = Binomial<T>(BinCoeff, i) - Binomial<T>(BinCoeff, i - 1);
            v(N - i - 1) = -v(i);
        }

        v /= std::pow(T(2), T(BinCoeff));
        return v;
    }
};

// Backward Hybrid Noise-Robust differentiators; http://www.holoborodko.com/pavel/wp-content/uploads/OneSidedNoiseRobustDifferentiators.pdf
template <typename T, int N> struct GetFHNRCoeffs {
    vectN_t<T, N> operator()() const;
};
template <typename T> struct GetFHNRCoeffs<T, 4> {
    vectN_t<T, 4> operator()() const { return (vectN_t<T, 4>() << T(2), T(-1), T(-2), T(1)).finished() / T(2); }
};
template <typename T> struct GetFHNRCoeffs<T, 5> {
    vectN_t<T, 5> operator()() const { return (vectN_t<T, 5>() << T(7), T(1), T(-10), T(-1), T(3)).finished() / T(10); }
};
template <typename T> struct GetFHNRCoeffs<T, 6> {
    vectN_t<T, 6> operator()() const { return (vectN_t<T, 6>() << T(16), T(1), T(-10), T(-10), T(-6), T(9)).finished() / T(28); }
};
template <typename T> struct GetFHNRCoeffs<T, 7> {
    vectN_t<T, 7> operator()() const { return (vectN_t<T, 7>() << T(12), T(5), T(-8), T(-6), T(-10), T(1), T(6)).finished() / T(28); }
};
template <typename T> struct GetFHNRCoeffs<T, 8> {
    vectN_t<T, 8> operator()() const { return (vectN_t<T, 8>() << T(22), T(7), T(-6), T(-11), T(-14), T(-9), T(-2), T(13)).finished() / T(60); }
};
template <typename T> struct GetFHNRCoeffs<T, 9> {
    vectN_t<T, 9> operator()() const { return (vectN_t<T, 9>() << T(52), T(29), T(-14), T(-17), T(-40), T(-23), T(-26), T(11), T(28)).finished() / T(180); }
};
template <typename T> struct GetFHNRCoeffs<T, 10> {
    vectN_t<T, 10> operator()() const { return (vectN_t<T, 10>() << T(56), T(26), T(-2), T(-17), T(-30), T(-30), T(-28), T(-13), T(4), T(34)).finished() / T(220); }
};
template <typename T> struct GetFHNRCoeffs<T, 11> {
    vectN_t<T, 11> operator()() const { return (vectN_t<T, 11>() << T(320), T(206), T(-8), T(-47), T(-186), T(-150), T(-214), T(-103), T(-92), T(94), T(180)).finished() / T(1540); }
};
template <typename T> struct GetFHNRCoeffs<T, 16> {
    vectN_t<T, 16> operator()() const { return (vectN_t<T, 16>() << T(322), T(217), T(110), T(35), T(-42), T(-87), T(-134), T(-149), T(-166), T(-151), T(-138), T(-93), T(-50), T(28), T(98), T(203)).finished() / T(2856); }
};

template <typename T, int N, typename BackwardCoeffs> vectN_t<T, N> GetBackwardISDCoeffs()
{
    vectN_t<T, N> v{};
    const vectN_t<T, N> v0 = BackwardCoeffs{}();
    for (Eigen::Index k = 0; k < N; ++k)
        v(k) = k * v0(k);
    return v;
}
// Backward Noise-Robust differentiators for irregular space data
template <typename T, int N> struct GetFNRISDCoeffs {
    vectN_t<T, N> operator()() const { return GetBackwardISDCoeffs<T, N, GetFNRCoeffs<T, N>>(); }
};
// Backward Hybrid Noise-Robust differentiators for irregular space data
template <typename T, int N> struct GetFHNRISDCoeffs {
    vectN_t<T, N> operator()() const { return GetBackwardISDCoeffs<T, N, GetFHNRCoeffs<T, N>>(); }
};

// Centered Noise-Robust differentiators (tangency at 2nd order): http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
template <typename T, int N>
struct GetCNR2Coeffs {
    vectN_t<T, N> operator()() const
    {
        static_assert(N % 2 == 1., "'N' must be odd.");
        return GetFNRCoeffs<T, N>{}(); // Same coefficients
    }
};

// Centered Noise-Robust  differentiators (tangency at 4th order): http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
template <typename T, int N> struct GetCNR4Coeffs {
    vectN_t<T, N> operator()() const;
};
template <typename T> struct GetCNR4Coeffs<T, 7> {
    vectN_t<T, 7> operator()() const { return (vectN_t<T, 7>() << T(-5), T(12), T(39), T(0), T(-39), T(-12), T(5)).finished() / T(96); }
};
template <typename T> struct GetCNR4Coeffs<T, 9> {
    vectN_t<T, 9> operator()() const { return (vectN_t<T, 9>() << T(-2), T(-1), T(16), T(27), T(0), T(-27), T(-16), T(1), T(2)).finished() / T(96); }
};
template <typename T> struct GetCNR4Coeffs<T, 11> {
    vectN_t<T, 11> operator()() const { return (vectN_t<T, 11>() << T(-11), T(-32), T(39), T(256), T(322), T(0), T(-322), T(-256), T(-39), T(32), T(11)).finished() / T(1536); }
};

// Centered Noise-Robust differentiators for irregular space data: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
template <typename T, int N, typename CNRCoeffs> vectN_t<T, N> GetCNRISDCoeffs()
{
    constexpr const int M = (N - 1) / 2;
    vectN_t<T, N> v{};
    const vectN_t<T, N> v0 = CNRCoeffs{}();
    v(M) = 0;
    for (int k = 1; k < M + 1; ++k) {
        v(M - k) = T(2) * k * v0(M - k);
        v(M + k) = T(2) * k * v0(M + k);
    }

    return v;
}
template <typename T, int N> struct GetCNR2ISDCoeffs {
    vectN_t<T, N> operator()() const { return GetCNRISDCoeffs<T, N, GetCNR2Coeffs<T, N>>(); }
};
template <typename T, int N> struct GetCNR4ISDCoeffs {
    vectN_t<T, N> operator()() const { return GetCNRISDCoeffs<T, N, GetCNR4Coeffs<T, N>>(); }
};

/*
 * Second order differentiators
 */

template <typename T>
constexpr T GetSONRBaseCoeff(int N, int M, int k)
{
    if (k > M)
        return T(0);
    else if (k == M)
        return T(1);
    else
        return ((T(2) * N - T(10)) * GetSONRBaseCoeff<T>(N, M, k + 1) - (N + T(2) * k + T(3)) * GetSONRBaseCoeff<T>(N, M, k + 2)) / (N - T(2) * k - T(1));
}

template <typename T, int N> vectN_t<T, (N - 1) / 2 + 1> GetSONRBaseCoeffs()
{
    static_assert(N >= 5 && N % 2 == 1, "N must be a odd number >= 5");
    constexpr const int M = (N - 1) / 2;
    vectN_t<T, M + 1> s{};
    for (int k = 0; k < M + 1; ++k)
        s(k) = GetSONRBaseCoeff<T>(N, M, k);

    return s;
}

// Second-Order Centered Noise-Robust differentiator: http://www.holoborodko.com/pavel/downloads/NoiseRobustSecondDerivative.pdf
template <typename T, int N>
struct GetSOCNRCoeffs {
    vectN_t<T, N> operator()() const
    {
        constexpr const int M = (N - 1) / 2;
        constexpr const T Den = pow(2, N - 3);
        vectN_t<T, N> v{};
        vectN_t<T, M + 1> s = GetSONRBaseCoeffs<T, N>();
        v(M) = s(0);
        for (int k = 1; k < M + 1; ++k) {
            v(M + k) = s(k);
            v(M - k) = s(k);
        }

        v /= Den;
        return v;
    }
};
// Second-Order Backward Noise-Robust differentiator: http://www.holoborodko.com/pavel/downloads/NoiseRobustSecondDerivative.pdf
template <typename T, int N> struct GetSOFNRCoeffs {
    vectN_t<T, N> operator()() const { return GetSOCNRCoeffs<T, N>(); }
}; // Coefficients are the same.

// Second-Order Centered Noise-Robust Irregular Space Data differentiator: http://www.holoborodko.com/pavel/downloads/NoiseRobustSecondDerivative.pdf
template <typename T, int N>
struct GetSOCNRISDCoeffs {
    vectN_t<T, N> operator()() const
    {
        constexpr const int M = (N - 1) / 2;
        constexpr const T Den = pow(2, N - 3);

        vectN_t<T, N> v{};
        const vectN_t<T, M + 1> s = GetSONRBaseCoeffs<T, N>();

        const auto alpha = [&s](int k) -> T { return T(4) * k * k * s(k); };

        v(M) = T(0);
        for (int k = 1; k < M + 1; ++k) {
            auto alph = alpha(k);
            v(M) -= T(2) * alph;
            v(M + k) = alph;
            v(M - k) = alph;
        }

        v /= Den;
        return v;
    }
};

// Second-Order Backward Noise-Robust Irregular Space Data differentiator: http://www.holoborodko.com/pavel/downloads/NoiseRobustSecondDerivative.pdf
template <typename T, int N> struct GetSOFNRISDCoeffs {
    vectN_t<T, N> operator()() const { return GetSOCNRISDCoeffs<T, N>(); }
}; // Same coefficients

/*
 * Differentiator Generator
 */

template <typename T, int N, int Order, typename CoeffGetter>
class BackwardDifferentiator : public GenericFilter<T> {
public:
    BackwardDifferentiator()
        : GenericFilter<T>(vectX_t<T>::Constant(1, T(1)), CoeffGetter{}())
    {}
    BackwardDifferentiator(T timestep)
        : GenericFilter<T>(vectX_t<T>::Constant(1, T(1)), CoeffGetter{}() / std::pow(timestep, Order))
    {}
    void setTimestep(T timestep) { this->setCoeffs(vectX_t<T>::Constant(1, T(1)), CoeffGetter{}() / std::pow(timestep, Order)); }
    T timestep() const noexcept { return std::pow(this->bCoeff()(0) / CoeffGetter{}()(0), T(1) / Order); }
};

template <typename T, int N, int Order, typename CoeffGetter>
class CenteredDifferentiator : public GenericFilter<T> {
public:
    CenteredDifferentiator()
        : GenericFilter<T>(vectX_t<T>::Constant(1, T(1)), CoeffGetter{}(), FilterType::Centered)
    {}
    CenteredDifferentiator(T timestep)
        : GenericFilter<T>(vectX_t<T>::Constant(1, T(1)), CoeffGetter{}() / std::pow(timestep, Order))
    {}
    void setTimestep(T timestep) { this->setCoeffs(vectX_t<T>::Constant(1, T(1)), CoeffGetter{}() / std::pow(timestep, Order)); }
    T timestep() const noexcept { return std::pow(this->bCoeff()(0) / CoeffGetter{}()(0), T(1) / Order); }
};

template <typename T, int N, int Order, typename CoeffGetter>
class TVBackwardDifferentiator : public TVGenericFilter<T> {
    static_assert(Order >= 1, "Order must be greater or equal to 1");

public:
    TVBackwardDifferentiator()
        : TVGenericFilter<T>(Order, vectX_t<T>::Constant(1, T(1)), CoeffGetter{}())
    {}
};

template <typename T, int N, int Order, typename CoeffGetter>
class TVCenteredDifferentiator : public TVGenericFilter<T> {
    static_assert(Order >= 1, "Order must be greater or equal to 1");

public:
    TVCenteredDifferentiator()
        : TVGenericFilter<T>(Order, vectX_t<T>::Constant(1, T(1)), CoeffGetter{}(), FilterType::Centered)
    {}
};

} // namespace details

// Backward differentiators
template <typename T, int N> using BackwardDiffNoiseRobust = details::BackwardDifferentiator<T, N, 1, details::GetFNRCoeffs<T, N>>;
template <typename T, int N> using BackwardDiffHybridNoiseRobust = details::BackwardDifferentiator<T, N, 1, details::GetFHNRCoeffs<T, N>>;
// Time-Varying backward differentiators
template <typename T, int N> using TVBackwardDiffNoiseRobust = details::TVBackwardDifferentiator<T, N, 1, details::GetFNRISDCoeffs<T, N>>;
template <typename T, int N> using TVBackwardDiffHybridNoiseRobust = details::TVBackwardDifferentiator<T, N, 1, details::GetFHNRISDCoeffs<T, N>>;

// Centered differentiators
template <typename T, int N> using CenteredDiffBasic = details::CenteredDifferentiator<T, N, 1, details::GetCDCoeffs<T, N>>;
template <typename T, int N> using CenteredDiffLowNoiseLanczos = details::CenteredDifferentiator<T, N, 1, details::GetLNLCoeffs<T, N>>;
template <typename T, int N> using CenteredDiffSuperLowNoiseLanczos = details::CenteredDifferentiator<T, N, 1, details::GetSLNLCoeffs<T, N>>;
template <typename T, int N> using CenteredDiffNoiseRobust2 = details::CenteredDifferentiator<T, N, 1, details::GetCNR2Coeffs<T, N>>;
template <typename T, int N> using CenteredDiffNoiseRobust4 = details::CenteredDifferentiator<T, N, 1, details::GetCNR4Coeffs<T, N>>;
// Time-Varying centered differentiators
template <typename T, int N> using TVCenteredDiffNoiseRobust2 = details::TVCenteredDifferentiator<T, N, 1, details::GetCNR2ISDCoeffs<T, N>>;
template <typename T, int N> using TVCenteredDiffNoiseRobust4 = details::TVCenteredDifferentiator<T, N, 1, details::GetCNR4ISDCoeffs<T, N>>;

// Second-order backward differentiators
template <typename T, int N> using BackwardDiffSecondOrder = details::BackwardDifferentiator<T, N, 2, details::GetSOFNRCoeffs<T, N>>;
// Second-order Time-Varying backward differentiators
template <typename T, int N> using TVBackwardDiffSecondOrder = details::TVBackwardDifferentiator<T, N, 2, details::GetSOFNRISDCoeffs<T, N>>;

// Second-order centered differentiators
template <typename T, int N> using CenteredDiffSecondOrder = details::CenteredDifferentiator<T, N, 2, details::GetSOCNRCoeffs<T, N>>;
// Second-order Time-Varying centered differentiators
template <typename T, int N> using TVCenteredDiffSecondOrder = details::TVCenteredDifferentiator<T, N, 2, details::GetSOCNRISDCoeffs<T, N>>;

} // namespace difi