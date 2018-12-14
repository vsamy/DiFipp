#include "BilinearTransform.h"
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
    updateCoeffSize();
    computeDigitalRep();
}

template <typename T>
void Butterworth<T>::computeDigitalRep()
{
    // Continuous pre-warped frequency
    T fpw = (m_fs / PI) * std::tan(PI * m_fc / m_fs);

    // Compute poles
    std::complex<T> analogPole;
    std::vector<std::complex<T>> poles(m_order);
    for (size_t k = 1; k <= m_order; ++k) {
        analogPole = generateAnalogPole(fpw, k);
        BilinearTransform<std::complex<T>>::SToZ(m_fs, analogPole, poles[k - 1]);
    }

    std::vector<std::complex<T>> zeros = generateAnalogZeros();
    std::vector<std::complex<T>> a = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(poles);
    std::vector<std::complex<T>> b = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(zeros);
    for (size_t i = 0; i < m_order + 1; ++i) {
        m_aCoeff[i] = a[i].real();
        m_bCoeff[i] = b[i].real();
    }

    scaleAmplitude();
    checkCoeff(m_aCoeff, m_bCoeff);
}

template <typename T>
void Butterworth<T>::updateCoeffSize()
{
    m_aCoeff.resize(m_order + 1);
    m_bCoeff.resize(m_order + 1);
    resetFilter();
}

template <typename T>
std::complex<T> Butterworth<T>::generateAnalogPole(T fpw, size_t k)
{
    T scaleFactor = 2 * PI * fpw;

    auto thetaK = [pi = PI, order = m_order](size_t k) -> T {
        return (2 * k - 1) * pi / (2 * order);
    };

    std::complex<T> analogPole(-std::sin(thetaK(k)), std::cos(thetaK(k)));
    switch (m_type) {
    case Type::HighPass:
        return scaleFactor / analogPole;

    case Type::LowPass:
    default:
        return scaleFactor * analogPole;
    }
}

template <typename T>
std::vector<std::complex<T>> Butterworth<T>::generateAnalogZeros()
{
    switch (m_type) {
    case Type::HighPass:
        return std::vector<std::complex<T>>(m_order, std::complex<T>(1));

    case Type::LowPass:
    default:
        return std::vector<std::complex<T>>(m_order, std::complex<T>(-1));
    }
}

template <typename T>
void Butterworth<T>::scaleAmplitude()
{
    T scale = 0;
    T sumB = 0;

    switch (m_type) {
    case Type::HighPass:
        for (size_t i = 0; i < m_order + 1; ++i) {
            if (i % 2 == 0) {
                scale += m_aCoeff[i];
                sumB += m_bCoeff[i];
            } else {
                scale -= m_aCoeff[i];
                sumB -= m_bCoeff[i];
            }
        }
        break;

    case Type::LowPass:
    default:
        for (size_t i = 0; i < m_order + 1; ++i) {
            scale += m_aCoeff[i];
            sumB += m_bCoeff[i];
        }
        break;
    }

    scale /= sumB;
    for (auto& b : m_bCoeff)
        b *= scale;
}

} // namespace fratio