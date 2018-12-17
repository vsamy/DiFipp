#include "BilinearTransform.h"
#include "polynome_functions.h"

namespace fratio {

template <typename T>
Butterworth<T>::Butterworth(Type type)
    : m_type(type)
{
}

template <typename T>
Butterworth<T>::Butterworth(int order, T fc, T fs, Type type)
    : m_type(type)
{
    initialize(order, fc, fs);
}

template <typename T>
void Butterworth<T>::setFilterParameters(int order, T fc, T fs)
{
    initialize(order, fc, fs);
}

template <typename T>
void Butterworth<T>::initialize(int order, T fc, T fs)
{
    if (order <= 0) {
        m_status = FilterStatus::BAD_ORDER_SIZE;
        return;
    }

    if (fc <= 0 || fs <= 0) {
        m_status = FilterStatus::BAD_FREQUENCY_VALUE;
        return;
    }

    if (m_fc > m_fs / 2.) {
        m_status = FilterStatus::BAD_CUTOFF_FREQUENCY;
        return;
    }

    m_order = order;
    m_fc = fc;
    m_fs = fs;
    computeDigitalRep();
    resetFilter();
}

template <typename T>
void Butterworth<T>::computeDigitalRep()
{
    // Continuous pre-warped frequency
    T fpw = (m_fs / PI) * std::tan(PI * m_fc / m_fs);

    // Compute poles
    std::complex<T> analogPole;
    vectX_t<std::complex<T>> poles(m_order);
    for (int k = 1; k <= m_order; ++k) {
        analogPole = generateAnalogPole(fpw, k);
        BilinearTransform<std::complex<T>>::SToZ(m_fs, analogPole, poles(k - 1));
    }

    vectX_t<std::complex<T>> zeros = generateAnalogZeros();
    vectX_t<std::complex<T>> a = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(poles);
    vectX_t<std::complex<T>> b = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(zeros);
    vectX_t<T> aCoeff(m_order + 1);
    vectX_t<T> bCoeff(m_order + 1);
    for (int i = 0; i < m_order + 1; ++i) {
        aCoeff(i) = a(i).real();
        bCoeff(i) = b(i).real();
    }

    scaleAmplitude(aCoeff, bCoeff);
    setCoeffs(std::move(aCoeff), std::move(bCoeff));
}

template <typename T>
std::complex<T> Butterworth<T>::generateAnalogPole(T fpw, int k)
{
    T scaleFactor = 2 * PI * fpw;

    auto thetaK = [pi = PI, order = m_order](int k) -> T {
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
vectX_t<std::complex<T>> Butterworth<T>::generateAnalogZeros()
{
    switch (m_type) {
    case Type::HighPass:
        return vectX_t<std::complex<T>>::Constant(m_order, std::complex<T>(1));

    case Type::LowPass:
    default:
        return vectX_t<std::complex<T>>::Constant(m_order, std::complex<T>(-1));
    }
}

template <typename T>
void Butterworth<T>::scaleAmplitude(Eigen::Ref<vectX_t<T>> aCoeff, Eigen::Ref<vectX_t<T>> bCoeff)
{
    T scale = 0;
    T sumB = 0;

    switch (m_type) {
    case Type::HighPass:
        for (int i = 0; i < m_order + 1; ++i) {
            if (i % 2 == 0) {
                scale += aCoeff(i);
                sumB += bCoeff(i);
            } else {
                scale -= aCoeff(i);
                sumB -= bCoeff(i);
            }
        }
        break;

    case Type::LowPass:
    default:
        scale = aCoeff.sum();
        sumB = bCoeff.sum();
        break;
    }

    bCoeff *= scale / sumB;
}

} // namespace fratio