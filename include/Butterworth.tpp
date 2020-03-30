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

#include "BilinearTransform.h"
#include "polynome_functions.h"

namespace difi {

template <typename T>
std::pair<int, T> Butterworth<T>::findMinimumButter(T wPass, T wStop, T APass, T AStop)
{
    Expects(wPass > T(0) && wPass < T(1));
    Expects(wStop > T(0) && wPass < T(1));
    T num = std::log10((std::pow(T(10), T(0.1) * std::abs(AStop)) - 1) / (std::pow(T(10), T(0.1) * std::abs(APass)) - 1));
    // pre-warp
    T fwPass = std::tan(T(0.5) * pi<T> * wPass);
    T fwStop = std::tan(T(0.5) * pi<T> * wStop);
    T w;
    if (wPass < wStop)
        w = std::abs(fwStop / fwPass);
    else
        w = std::abs(fwPass / fwStop);
    T denum = T(2) * std::log10(w);

    int order = static_cast<int>(std::ceil(num / denum));
    T ctf = w / std::pow(std::pow(T(10), T(0.1) * std::abs(AStop)) - 1, T(1) / T(2 * order));
    if (wPass < wStop)
        ctf *= fwPass;
    else
        ctf = fwPass / ctf;

    return std::pair<int, T>(order, T(2) * std::atan(ctf) / pi<T>);
}

template <typename T>
Butterworth<T>::Butterworth(Type type)
    : m_type(type)
{
}

template <typename T>
Butterworth<T>::Butterworth(int order, T fc, T fs, Type type)
    : m_type(type)
{
    setFilterParameters(order, fc, fs);
}

template <typename T>
Butterworth<T>::Butterworth(int order, T fLower, T fUpper, T fs, Type type)
    : m_type(type)
{
    setFilterParameters(order, fLower, fUpper, fs);
}

template <typename T>
void Butterworth<T>::setFilterParameters(int order, T fc, T fs)
{
    Expects(fc < fs / T(2));
    initialize(order, fc, 0, fs);
}

template <typename T>
void Butterworth<T>::setFilterParameters(int order, T fLower, T fUpper, T fs)
{
    Expects(fLower < fUpper);
    initialize(order, fLower, fUpper, fs);
}

template <typename T>
void Butterworth<T>::initialize(int order, T f1, T f2, T fs)
{
    // f1 = fc for LowPass/HighPass filter
    // f1 = fLower, f2 = fUpper for BandPass/BandReject filter
    Expects(order > 0);
    Expects(f1 > 0 && fs > 0); // f2 must be > f1 check in setFilterParameters

    m_order = order;
    m_fs = fs;
    if (m_type == Type::LowPass || m_type == Type::HighPass)
        computeDigitalRep(f1);
    else
        computeBandDigitalRep(f1, f2); // For band-like filters
}

template <typename T>
void Butterworth<T>::computeDigitalRep(T fc)
{
    // Continuous pre-warped frequency
    T fpw = (m_fs / pi<T>)*std::tan(pi<T> * fc / m_fs);

    // Compute poles
    vectXc_t<T> poles(m_order);
    std::complex<T> analogPole;
    for (int k = 0; k < m_order; ++k) {
        analogPole = generateAnalogPole(k + 1, fpw);
        BilinearTransform<std::complex<T>>::SToZ(m_fs, analogPole, poles(k));
    }

    vectXc_t<T> zeros = generateAnalogZeros();
    vectXc_t<T> a = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(poles);
    vectXc_t<T> b = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(zeros);
    vectX_t<T> aCoeff(m_order + 1);
    vectX_t<T> bCoeff(m_order + 1);
    for (int i = 0; i < m_order + 1; ++i) {
        aCoeff(i) = a(i).real();
        bCoeff(i) = b(i).real();
    }

    scaleAmplitude(aCoeff, bCoeff);
    this->setCoeffs(std::move(aCoeff), std::move(bCoeff));
}

template <typename T>
void Butterworth<T>::computeBandDigitalRep(T fLower, T fUpper)
{
    T fpw1 = (m_fs / pi<T>)*std::tan(pi<T> * fLower / m_fs);
    T fpw2 = (m_fs / pi<T>)*std::tan(pi<T> * fUpper / m_fs);
    T fpw0 = std::sqrt(fpw1 * fpw2);

    vectXc_t<T> poles(2 * m_order);
    std::pair<std::complex<T>, std::complex<T>> analogPoles;
    for (int k = 0; k < m_order; ++k) {
        analogPoles = generateBandAnalogPole(k + 1, fpw0, fpw2 - fpw1);
        BilinearTransform<std::complex<T>>::SToZ(m_fs, analogPoles.first, poles(k));
        BilinearTransform<std::complex<T>>::SToZ(m_fs, analogPoles.second, poles(m_order + k));
    }

    vectXc_t<T> zeros = generateAnalogZeros(fpw0);
    vectXc_t<T> a = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(poles);
    vectXc_t<T> b = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(zeros);
    vectX_t<T> aCoeff(2 * m_order + 1);
    vectX_t<T> bCoeff(2 * m_order + 1);
    for (int i = 0; i < 2 * m_order + 1; ++i) {
        aCoeff(i) = a(i).real();
        bCoeff(i) = b(i).real();
    }

    if (m_type == Type::BandPass)
        scaleAmplitude(aCoeff, bCoeff, std::exp(std::complex<T>(T(0), T(2) * pi<T> * std::sqrt(fLower * fUpper) / m_fs)));
    else
        scaleAmplitude(aCoeff, bCoeff);

    this->setCoeffs(std::move(aCoeff), std::move(bCoeff));
}

template <typename T>
std::complex<T> Butterworth<T>::generateAnalogPole(int k, T fpw1)
{
    auto thetaK = [pi = pi<T>, order = m_order](int k) -> T {
        return static_cast<float>(2 * k - 1) * pi / static_cast<float>(2 * order);
    };

    std::complex<T> analogPole(-std::sin(thetaK(k)), std::cos(thetaK(k)));
    switch (m_type) {
    case Type::HighPass:
        return T(2) * pi<T> * fpw1 / analogPole;
    case Type::LowPass:
        return T(2) * pi<T> * fpw1 * analogPole;
    default:
        GSL_ASSUME(0);
    }
}

template <typename T>
std::pair<std::complex<T>, std::complex<T>> Butterworth<T>::generateBandAnalogPole(int k, T fpw0, T bw)
{
    auto thetaK = [pi = pi<T>, order = m_order](int k) -> T {
        return static_cast<float>(2 * k - 1) * pi / static_cast<float>(2 * order);
    };

    std::complex<T> analogPole(-std::sin(thetaK(k)), std::cos(thetaK(k)));
    std::pair<std::complex<T>, std::complex<T>> poles;
    std::complex<T> s0 = T(2) * pi<T> * fpw0;
    std::complex<T> s = T(0.5) * bw / fpw0;
    switch (m_type) {
    case Type::BandReject:
        s /= analogPole;
        poles.first = s0 * (s + std::complex<T>(T(0), T(1)) * std::sqrt(T(1) - s * s));
        poles.second = s0 * (s - std::complex<T>(T(0), T(1)) * std::sqrt(T(1) - s * s));
        return poles;
    case Type::BandPass:
        s *= analogPole;
        poles.first = s0 * (s + std::complex<T>(T(0), T(1)) * std::sqrt(T(1) - s * s));
        poles.second = s0 * (s - std::complex<T>(T(0), T(1)) * std::sqrt(T(1) - s * s));
        return poles;
    default:
        GSL_ASSUME(0);
    }
}

template <typename T>
vectXc_t<T> Butterworth<T>::generateAnalogZeros(T fpw0)
{
    switch (m_type) {
    case Type::HighPass:
        return vectXc_t<T>::Constant(m_order, std::complex<T>(1));
    case Type::BandPass:
        return (vectXc_t<T>(2 * m_order) << vectXc_t<T>::Constant(m_order, std::complex<T>(-1)), vectXc_t<T>::Constant(m_order, std::complex<T>(1))).finished();
    case Type::BandReject: {
        T w0 = T(2) * std::atan(pi<T> * fpw0 / m_fs);
        return (vectXc_t<T>(2 * m_order) << vectXc_t<T>::Constant(m_order, std::exp(std::complex<T>(0, w0))), vectXc_t<T>::Constant(m_order, std::exp(std::complex<T>(0, -w0)))).finished();
    }
    case Type::LowPass:
    default:
        return vectXc_t<T>::Constant(m_order, std::complex<T>(-1));
    }
}

template <typename T>
void Butterworth<T>::scaleAmplitude(const vectX_t<T>& aCoeff, Eigen::Ref<vectX_t<T>> bCoeff, const std::complex<T>& bpS)
{
    T num = 0;
    T denum = 0;

    switch (m_type) {
    case Type::HighPass:
        for (int i = 0; i < m_order + 1; ++i) {
            if (i % 2 == 0) {
                num += aCoeff(i);
                denum += bCoeff(i);
            } else {
                num -= aCoeff(i);
                denum -= bCoeff(i);
            }
        }
        break;
    case Type::BandPass: {
        std::complex<T> numC(bCoeff(0));
        std::complex<T> denumC(aCoeff(0));
        for (int i = 1; i < 2 * m_order + 1; ++i) {
            denumC = denumC * bpS + aCoeff(i);
            numC = numC * bpS + bCoeff(i);
        }
        num = std::abs(denumC);
        denum = std::abs(numC);
    } break;
    case Type::BandReject:
    case Type::LowPass:
    default:
        num = aCoeff.sum();
        denum = bCoeff.sum();
        break;
    }

    bCoeff *= num / denum;
}

} // namespace difi