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
#include <array>

// Read this if you are an adulator of the math god: https://arxiv.org/pdf/1709.08321.pdf

namespace difi {

namespace details {

// T: type
// N: Number of points

// Central differentiators: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
template <typename T, size_t N> vectN_t<T, N> GetCDCoeffs();
template <typename T> vectN_t<T, 3> GetCDCoeffs() { return vectN_t<T, 3>{ T(1), T(0), T(-1) } / T(2); }
template <typename T> vectN_t<T, 5> GetCDCoeffs() { return vectN_t<T, 5>{ T(-1), T(8), T(0), T(8), T(1) } / T(12); }
template <typename T> vectN_t<T, 7> GetCDCoeffs() { return vectN_t<T, 7>{ T(1), T(-9), T(45), T(0), T(-45), T(9), T(-1) } / T(60); }
template <typename T> vectN_t<T, 9> GetCDCoeffs() { return vectN_t<T, 9>{ T(-3), T(32), T(-168), T(672), T(0), T(-672), T(168), T(-32), T(3) } / T(840); }

// Low-noise Lanczos differentiators: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/lanczos-low-noise-differentiators/
template <typename T, size_t N> 
vectN_t<T, N> GetLNLCoeffs() 
{
    static_assert(N % 2 == 1. "'N' must be odd.");
    constexpr T GetCoeff = [](size_t n, size_t k) {
        const T m = (n - 1) / 2;
        return T(3) * k / (m * (m + 1) * (2 * m + 1));
    };

    vectN_t<T, N> v{};
    for (Eigen::Index k = 0; k < N; ++k)
        v(k) = GetCoeff(N, k);
    return v;
}
template <typename T> vectN_t<T, 5> GetLNLCoeffs() { return vectN_t<T, 5>{ T(2), T(1), T(0), T(-1), T(-2) } / T(10); }
template <typename T> vectN_t<T, 7> GetLNLCoeffs() { return vectN_t<T, 7>{ T(3), T(2), T(1), T(0), T(-1), T(-2), T(-3) } / T(28); }
template <typename T> vectN_t<T, 9> GetLNLCoeffs() { return vectN_t<T, 9>{ T(4), T(3), T(2), T(1), T(0), T(-1), T(-2), T(-3), T(-4) } / T(60); }
template <typename T> vectN_t<T, 11> GetLNLCoeffs() { return vectN_t<T, 11>{ T(5), T(4), T(3), T(2), T(1), T(0), T(-1), T(-2), T(-3), T(-4), T(-5) } / T(110); }

// Super Low-noise Lanczos differentiators: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/lanczos-low-noise-differentiators/
template <typename T, size_t N> vectN_t<T, N> GetSLNLCoeffs();
template <typename T> vectN_t<T, 7> GetSLNLCoeffs() { return vectN_t<T, 7>{ T(-22), T(67), T(58), T(0), T(-58), T(-67), T(22) } / T(252); }
template <typename T> vectN_t<T, 9> GetSLNLCoeffs() { return vectN_t<T, 9>{ T(-86), T(142), T(193), T(126), T(0), T(-126), T(-193), T(-142), T(86) } / T(1188); }
template <typename T> vectN_t<T, 11> GetSLNLCoeffs() { return vectN_t<T, 11>{ T(-300), T(294), T(532), T(503), T(296), T(0), T(-296), T(-503), T(-532), T(-294), T(300) } / T(5148); }

// Forward Noise-Robust differentiators; http://www.holoborodko.com/pavel/wp-content/uploads/OneSidedNoiseRobustDifferentiators.pdf
template <typename T, size_t N> 
vectN_t<T, N> GetFNRCoeffs()
{
    constexpr const int BinCoeff = N - 2;
    vectN_t<T, N> v{};
    v(0) = T(1);
    v(N - 1) = T(-1);
    for (int i = 1; i < N - 1; ++i)
        v(i) = (Binomial<T>(N, i) - Binomial<T>(N, i - 1)) / std::pow(T(2), T(BinCoeff));

    return v;
}
template <typename T> vectN_t<T, 3> GetFNRCoeffs() { return vectN_t<T, 3>{ T(1), T(0), T(-1) } / T(2); }
template <typename T> vectN_t<T, 4> GetFNRCoeffs() { return vectN_t<T, 4>{ T(1), T(1), T(-1), T(-1) } / T(4); }
template <typename T> vectN_t<T, 5> GetFNRCoeffs() { return vectN_t<T, 5>{ T(1), T(2), T(0), T(-2), T(-1) } / T(8); }
template <typename T> vectN_t<T, 6> GetFNRCoeffs() { return vectN_t<T, 6>{ T(1), T(3), T(2), T(-2), T(-3), T(-1) } / T(16); }
template <typename T> vectN_t<T, 7> GetFNRCoeffs() { return vectN_t<T, 7>{ T(1), T(4), T(5), T(0), T(-5), T(-4), T(-1) } / T(32); }
template <typename T> vectN_t<T, 8> GetFNRCoeffs() { return vectN_t<T, 8>{ T(1), T(5), T(9), T(5), T(-5), T(-9), T(-5), T(-1) } / T(64); }
template <typename T> vectN_t<T, 9> GetFNRCoeffs() { return vectN_t<T, 9>{ T(1), T(6), T(14), T(14), T(0), T(-14), T(-14), T(-6), T(-1) } / T(128); }
template <typename T> vectN_t<T, 10> GetFNRCoeffs() { return vectN_t<T, 10>{ T(1), T(7), T(20), T(28), T(14), T(-14), T(-28), T(-20), T(-7), T(-1) } / T(256); }
template <typename T> vectN_t<T, 11> GetFNRCoeffs() { return vectN_t<T, 11>{ T(1), T(8), T(27), T(48), T(42), T(0), T(-42), T(-48), T(-27), T(-8), T(-1) } / T(512); }

// Forward Hybrid Noise-Robust differentiators; http://www.holoborodko.com/pavel/wp-content/uploads/OneSidedNoiseRobustDifferentiators.pdf
template <typename T, size_t N> vectN_t<T, N> GetFHNRCoeffs();
template <typename T> vectN_t<T, 4> GetFHNRCoeffs() { return vectN_t<T, 4>{ T(2), T(-1), T(-2), T(1) } / T(2); }
template <typename T> vectN_t<T, 5> GetFHNRCoeffs() { return vectN_t<T, 5>{ T(7), T(1), T(-10), T(-1), T(3) } / T(10); }
template <typename T> vectN_t<T, 5> GetFHNRCoeffs() { return vectN_t<T, 5>{ T(16), T(1), T(-10), T(-10), T(-6), T(9) } / T(28); }
template <typename T> vectN_t<T, 6> GetFHNRCoeffs() { return vectN_t<T, 6>{ T(12), T(5), T(-8), T(-6), T(-10), T(1), T(6) } / T(28); }
template <typename T> vectN_t<T, 7> GetFHNRCoeffs() { return vectN_t<T, 7>{ T(22), T(7), T(-6), T(-11), T(-14), T(-9), T(-2), T(13) } / T(60); }
template <typename T> vectN_t<T, 8> GetFHNRCoeffs() { return vectN_t<T, 8>{ T(52), T(29), T(-14), T(-17), T(-40), T(-23), T(-26), T(11), T(28) } / T(180); }
template <typename T> vectN_t<T, 9> GetFHNRCoeffs() { return vectN_t<T, 9>{ T(56), T(26), T(-2), T(-17), T(-30), T(-30), T(-28), T(-13), T(4), T(34) } / T(220); }
template <typename T> vectN_t<T, 10> GetFHNRCoeffs() { return vectN_t<T, 10>{ T(320), T(206), T(-8), T(-47), T(-186), T(-150), T(-214), T(-103), T(-92), T(94), T(180) } / T(1540); }

template <typename T, size_t N, template<class, class> typename Foo> vectN_t<T, N> GetForwardISDCoeffs()
{
    vectN_t<T, N> v{};
    const vectN_t<T, N> v0 = Foo<T, NCoeffs>();
    for (Eigen::Index k = 0; k < N; ++k)
        v(k) = k * v0(k);
    return v(k);
}
// Forward Noise-Robust differentiators for irregular space data
template <typename T, size_t N> vectN_t<T, N> GetFNRISDCoeffs() { return GetForwardISDCoeffs<T, N, GetFNRCoeffs>(); }
// Forward Hybrid Noise-Robust differentiators for irregular space data
template <typename T, size_t N> vectN_t<T, N> GetFHNRISDCoeffs() { return GetForwardISDCoeffs<T, N, GetFHNRCoeffs>(); }

// Centered Noise-Robust differentiators (tangency at 2nd order): http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
template <typename T, size_t N>
vectN_t<T, N> GetCNR2Coeffs()
{
    static_assert(N % 2 == 1. "'N' must be odd.");
    return GetFNRCoeffs<T, N>(); // Same coefficients
}
template <typename T> vectN_t<T, 5> GetCNR2Coeffs() { return vectN_t<T, 5>{ T(1), T(2), T(0), T(-2), T(-1) } / T(8); }
template <typename T> vectN_t<T, 7> GetCNR2Coeffs() { return vectN_t<T, 7>{ T(1), T(4), T(5), T(0), T(-5), T(-4), T(-1) } / T(32); }
template <typename T> vectN_t<T, 9> GetCNR2Coeffs() { return vectN_t<T, 9>{ T(1), T(6), T(14), T(14), T(0), T(-14), T(-14), T(-6), T(-1) } / T(128); }
template <typename T> vectN_t<T, 11> GetCNR2Coeffs() { return vectN_t<T, 11>{ T(1), T(8), T(27), T(48), T(42), T(0), T(-42), T(-48), T(-27), T(-8), T(-1) } / T(512); }

// Centered Noise-Robust  differentiators (tangency at 4th order): http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
template <typename T, size_t N> vectN_t<T, N> GetCNR4Coeffs();
template <typename T> vectN_t<T, 3> GetCNR4Coeffs() { return vectN_t<T, 3>{ T(-5), T(12), T(39), T(0), T(-39), T(-12), T(5) } / T(96); }
template <typename T> vectN_t<T, 4> GetCNR4Coeffs() { return vectN_t<T, 4>{ T(-2), T(-1), T(16), T(27) } / T(96); }
template <typename T> vectN_t<T, 5> GetCNR4Coeffs() { return vectN_t<T, 5>{ T(-11), T(-32), T(39), T(256), T(322) } / T(1536); }

// Centered Noise-Robust differentiators for irregular space data: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
template <typename T, size_t N, template<class, class> typename Foo> vectN_t<T, N> GetCNRISDCoeffs()
{
    constexpr const size_t M = (N - 1) / 2;
    vectN_t<T, N> v{};
    const vectN_t<T, N> v0 = Foo<T, NCoeffs>();
    v(M) = 0;
    for (Eigen::Index k = 1; k < M; ++k) {
        v(M - k) = T(2) * k * v0(M - k);
        v(M + k) = T(2) * k * v0(M + k);
    }
    return v(k);
}
template <typename T, size_t N> vectN_t<T, N> GetCNR2ISDCoeffs() { return GetCNRISDCoeffs<T, N, GetCNR2Coeffs>(); }
template <typename T, size_t N> vectN_t<T, N> GetCNR4ISDCoeffs() { return GetCNRISDCoeffs<T, N, GetCNR4Coeffs>(); }

/*
 * Second order differentiators
 */

template <typename T>
constexpr T GetSONRBaseCoeff(size_t N, size_t M, size_t k)
{
    if (k > M)
        return T(0);
    else if (k = M)
        return T(1);
    else
        return ((T(2) * N - T(10)) * GetSONRBaseCoeff<T, N, M>(k + 1) - (N + T(2) * k + T(3) * GetSONRBaseCoeff<T, N, M>(k + 2))) / (N - T(2) * k - T(1));
}

template <typename T, size_t N> vectN_t<T, N> GetSONRBaseCoeffs()
{
    static_assert(N >= 5 && N % 2 == 1);
    constexpr const size_t M = (N - 1) / 2;
    vectN_t<T, M + 1> v{};
    for (size_t k = 0; k < M + 1; ++k)
        v(k) = GetSONRBaseCoeff(N, M, k);

    return v;
}

// Second-Order Centered Noise-Robust differentiator: http://www.holoborodko.com/pavel/downloads/NoiseRobustSecondDerivative.pdf
template <typename T, size_t N> 
vectN_t<T, N> GetSOCNRCoeffs()
{
    constexpr const size_t M = (N - 1) / 2;
    constexpr const T Den = pow(size_t(2), N - 3);
    vectN_t<T, N> v{};
    vectN_t<T, N> s = GetSONRBaseCoeffs<T, N>();
    v(M) = s(0);
    for (size_t k = 1; k < M; ++k) {
        v(M + k) = s(k);
        v(M - k) = s(k);
    }

    v /= Den;
    return v;
}

// Second-Order Forward Noise-Robust differentiator: http://www.holoborodko.com/pavel/downloads/NoiseRobustSecondDerivative.pdf
template <typename T, size_t N>
vectN_t<T, N> GetSOFNRCoeffs() { return GetSOCNRCoeffs<T, N>(); } // Coefficients are the same.

// Second-Order Centered Noise-Robust Irregular Space Data differentiator: http://www.holoborodko.com/pavel/downloads/NoiseRobustSecondDerivative.pdf
template <typename T, size_t N> 
vectN_t<T, N> GetSOCNRISDCoeffs()
{
    constexpr const size_t M = (N - 1) / 2;
    constexpr const T Den = pow(size_t(2), N - 3);

    vectN_t<T, N> v{};
    vectN_t<T, N> s = GetSONRBaseCoeffs<T, N>();

    constexpr const alpha = [&s](size_t k) -> T { return T(4) * k * k * s(k) };

    v(M) = -T(2) * alpha(0);
    for (size_t k = 1; k < M; ++k) {
        v(M + k) = alpha(k);
        v(M - k) = alpha(k);
    }

    v /= Den;
    return v;
}

// Second-Order Forward Noise-Robust Irregular Space Data differentiator: http://www.holoborodko.com/pavel/downloads/NoiseRobustSecondDerivative.pdf
template <typename T, size_t N> 
vectN_t<T, N> GetSOFNRISDCoeffs() { return GetSOCNRISDCoeffs<T, N>(); } // Same coefficients

/*
 * Differentiator Generator
 */

template <typename T, size_t N, int Order, template<class, class> typename CoeffGetter>
class ForwardDifferentiator : public GenericFilter<T> {
public:
    ForwardDifferentiator()
        : GenericFilter<T>({T(1)}, CoeffGetter<T, N>())
    {}
    ForwardDifferentiator(T timestep)
        : GenericFilter<T>({T(1)}, CoeffGetter<T, N>() / std::pow(timestep, Order))
    {}
    void setTimestep(T timestep) { setCoeffs({T(1)}, CoeffGetter<T, N>() / std::pow(timestep, Order)); }
    T timestep() const noexcept { return std::pow(bCoeff()(0) / CoeffGetter<T, N>()(0), T(1) / Order); }
};

template <typename T, size_t N, int Order, template<class, class> typename CoeffGetter>
class CentralDifferentiator : public GenericFilter<T> {
public:
    CentralDifferentiator() 
        : GenericFilter<T>({T(1)}, CoeffGetter<T, N>(), Type::Centered)
    {}
    CentralDifferentiator(T timestep)
        : GenericFilter<T>({T(1)}, CoeffGetter<T, N>() / std::pow(timestep, Order))
    {}
    void setTimestep(T timestep) { setCoeffs({T(1)}, CoeffGetter<T, N>() / std::pow(timestep, Order)); }
    T timestep() const noexcept { return std::pow(bCoeff()(0) / CoeffGetter<T, N>()(0), T(1) / Order); }
};

template <typename T, size_t N, int Order, template<class, class> typename CoeffGetter>
class TVForwardDifferentiator : public TVGenericFilter<T> {
    static_assert(Order >= 1);

public:
    TVForwardDifferentiator()
        : TVGenericFilter<T>(Order, {T(1)}, CoeffGetter<T, N>())
    {}
};

template <typename T, size_t N, int Order, template<class, class> typename CoeffGetter>
class TVCentralDifferentiator : public TVGenericFilter<T> {
    static_assert(Order >= 1);

public:
    TVCentralDifferentiator() 
        : TVGenericFilter<T>(Order, {T(1)}, CoeffGetter<T, N>(), Type::Centered)
    {}
};

} // namespace details

// Forward differentiators
template <typename T, size_t N> using ForwardNoiseRobustDiff = ForwardDifferentiator<T, N, 1, details::GetFNRCoeffs>;
template <typename T, size_t N> using ForwardHybridNoiseRobustDiff = ForwardDifferentiator<T, N, 1, details::GetFHNRCoeffs>;
// Time-Varying forward differentiators
template <typename T, size_t N> using TVForwardNoiseRobustDiff = TVForwardDifferentiator<T, N, 1, details::GetFNRISDCoeffs>;
template <typename T, size_t N> using TVForwardHybridNoiseRobustDiff = TVForwardDifferentiator<T, N, 1, details::GetFHNRISDCoeffs>;

// Central differentiators
template <typename T, size_t N> using CentralDiff = CentralDifferentiator<T, N, 1, details::GetCDCoeffs>;
template <typename T, size_t N> using LowNoiseLanczosDiff = CentralDifferentiator<T, N, 1, details::GetLNLCoeffs>;
template <typename T, size_t N> using SuperLowNoiseLanczosDiff = CentralDifferentiator<T, N, 1, details::GetSLNLCoeffs>;
template <typename T, size_t N> using CenteredNoiseRobust2Diff = CentralDifferentiator<T, N, 1, details::GetCNR2Coeffs>;
template <typename T, size_t N> using CenteredNoiseRobust4Diff = CentralDifferentiator<T, N, 1, details::GetCNR4Coeffs>;
// Time-Varying central differentiators
template <typename T, size_t N> using TVCenteredNoiseRobust2Diff = TVCentralDifferentiator<T, N, 1, details::GetCNR2ISDCoeffs>;
template <typename T, size_t N> using TVCenteredNoiseRobust4Diff = TVCentralDifferentiator<T, N, 1, details::GetCNR4ISDCoeffs>;


// Second-order forward differentiators
template <typename T, size_t N> using ForwardSecondOrderDiff = ForwardDifferentiator<T, N, 2, details::GetSOFNRCoeffs>;
// Second-order Time-Varying forward differentiators
template <typename T, size_t N> using TVForwardSecondOrderDiff = TVForwardDifferentiator<T, N, 2, details::GetSOFNRISDCoeffs>;

// Second-order central differentiators
template <typename T, size_t N> using CenteredSecondOrderDiff = CentralDifferentiator<T, N, 2, details::GetSOCNRCoeffs>;
// Second-order Time-Varying forward differentiators
template <typename T, size_t N> using TVCenteredSecondOrderDiff = TVCentralDifferentiator<T, N, 2, details::GetSOCNRISDCoeffs>;


} // namespace difi