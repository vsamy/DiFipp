#include "polynome_functions.h"
#include <cmath>
#include <sstream>

namespace fratio {

template <typename T>
Butterworth<T>::Butterworth(Type type)
    : m_type(type)
{
}

template <typename T>
Butterworth<T>::Butterworth(size_t order, T fc, T fs, Type type)
    : m_type(type)
{
    initialize(order, fc, fs);
}

template <typename T>
void Butterworth<T>::setFilterParameters(size_t order, T fc, T fs)
{
    initialize(order, fc, fs);
}

template <typename T>
void Butterworth<T>::initialize(size_t order, T fc, T fs)
{
    if (m_fc > m_fs / 2.) {
        std::stringstream ss;
        ss << "The cut-off-frequency must be inferior to the sampling frequency"
           << "\n  Given cut-off-frequency is " << m_fc
           << "\n  Given sample frequency is " << m_fs;
        throw std::runtime_error(ss.str());
    }

    m_order = order;
    m_fc = fc;
    m_fs = fs;
    m_poles.resize(order);
    updateCoeffSize();
    computeDigitalRep();
}

template <typename T>
void Butterworth<T>::computeDigitalRep()
{
    T pi = static_cast<T>(M_PI);
    // Continuous pre-warped frequency
    T fpw = (m_fs / pi) * std::tan(pi * m_fc / m_fs);
    T scaleFactor = T(2) * pi * fpw;

    auto thetaK = [pi, order = m_order](size_t k) -> T {
        return (T(2) * k - T(1)) * pi / (T(2) * order);
    };

    // Compute poles
    std::complex<T> scalePole;
    for (size_t k = 1; k <= m_order; ++k) {
        scalePole = scaleFactor * std::complex<T>(-std::sin(thetaK(k)), std::cos(thetaK(k)));
        scalePole /= T(2) * m_fs;
        m_poles[k - 1] = (T(1) + scalePole) / (T(1) - scalePole);
    }

    std::vector<std::complex<T>> numPoles(m_order, std::complex<T>(-1));
    std::vector<std::complex<T>> a = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(m_poles);
    std::vector<std::complex<T>> b = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(numPoles);
    T norm = 0;
    T sumB = 0;
    for (size_t i = 0; i < m_order + 1; ++i) {
        m_aCoeff[i] = a[i].real();
        m_bCoeff[i] = b[i].real();
        norm += m_aCoeff[i];
        sumB += m_bCoeff[i];
    }

    norm /= sumB;
    for (auto& b : m_bCoeff)
        b *= norm;

    checkCoeff(m_aCoeff, m_bCoeff);
}

template <typename T>
void Butterworth<T>::updateCoeffSize()
{
    m_aCoeff.resize(m_order + 1);
    m_bCoeff.resize(m_order + 1);
    resetFilter();
}

} // namespace fratio