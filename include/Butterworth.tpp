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
Butterworth<T>::Butterworth(int order, T fc, T fs, T fCenter, Type type)
    : m_type(type)
{
    initialize(order, fc, fs, fCenter);
}

template <typename T>
void Butterworth<T>::setFilterParameters(int order, T fc, T fs, T fCenter)
{
    initialize(order, fc, fs, fCenter);
}

template <typename T>
void Butterworth<T>::initialize(int order, T fc, T fs, T fCenter)
{
    if (order <= 0) {
        m_status = FilterStatus::BAD_ORDER_SIZE;
        return;
    }

    if (fc <= 0 || fs <= 0) {
        m_status = FilterStatus::BAD_FREQUENCY_VALUE;
        return;
    }

    if ((m_type == Type::BandPass || m_type == Type::BandReject) && fCenter - fc / 2. <= 0) {
        m_status = FilterStatus::BAD_BAND_FREQUENCY;
        return;
    }

    if (m_fc > m_fs / 2.) {
        m_status = FilterStatus::BAD_CUTOFF_FREQUENCY;
        return;
    }

    m_order = order;
    m_fs = fs;
    if (m_type == Type::LowPass || m_type == Type::HighPass)
        computeDigitalRep(fc);
    else
        computeBandDigitalRep(fc, fCenter); // For band-like filters

    resetFilter();
}

template <typename T>
void Butterworth<T>::computeDigitalRep(T fc)
{
    m_fc = fc;
    // Continuous pre-warped frequency
    T fpw = (m_fs / PI) * std::tan(PI * m_fc / m_fs);

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
    setCoeffs(std::move(aCoeff), std::move(bCoeff));
}

template <typename T>
void Butterworth<T>::computeBandDigitalRep(T bw, T fCenter)
{
    m_bw = bw;
    m_fCenter = fCenter;
    T f1 = m_fCenter - m_bw / T(2);
    T f2 = m_fCenter + m_bw / T(2);
    T fpw1 = (m_fs / PI) * std::tan(PI * f1 / m_fs);
    T fpw2 = (m_fs / PI) * std::tan(PI * f2 / m_fs);

    vectXc_t<T> poles(2 * m_order);
    std::pair<std::complex<T>, std::complex<T>> analogPoles;
    for (int k = 0; k < m_order; ++k) {
        analogPoles = generateBandAnalogPole(k + 1, fpw1, fpw2);
        BilinearTransform<std::complex<T>>::SToZ(m_fs, analogPoles.first, poles(2 * k));
        BilinearTransform<std::complex<T>>::SToZ(m_fs, analogPoles.second, poles(2 * k + 1));
    }

    vectXc_t<T> zeros = generateAnalogZeros();
    vectXc_t<T> a = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(poles);
    vectXc_t<T> b = VietaAlgo<std::complex<T>>::polyCoeffFromRoot(zeros);
    vectX_t<T> aCoeff(2 * m_order + 1);
    vectX_t<T> bCoeff(2 * m_order + 1);
    for (int i = 0; i < 2 * m_order + 1; ++i) {
        aCoeff(i) = a(i).real();
        bCoeff(i) = b(i).real();
    }

    std::complex<T> s = std::exp(std::complex<T>(T(0), T(2) * PI * std::sqrt(f1 * f2) / m_fs));
    std::complex<T> num(b(0));
    std::complex<T> denum(a(0));
    for (int i = 1; i < 2 * m_order + 1; ++i) {
        num = num * s + b(i);
        denum = denum * s + a(i);
    }
    std::complex<T> hf0 = num / denum;
    bCoeff *= T(1) / std::abs(hf0); // bCoeff *= 1 / abs(num / denum)
    setCoeffs(std::move(aCoeff), std::move(bCoeff));
}

template <typename T>
std::complex<T> Butterworth<T>::generateAnalogPole(int k, T fpw1)
{
    auto thetaK = [pi = PI, order = m_order](int k) -> T {
        return (2 * k - 1) * pi / (2 * order);
    };

    std::complex<T> analogPole(-std::sin(thetaK(k)), std::cos(thetaK(k)));
    switch (m_type) {
    case Type::HighPass:
        return T(2) * PI * fpw1 / analogPole;
    case Type::LowPass:
    default:
        return T(2) * PI * fpw1 * analogPole;
    }
}

template <typename T>
std::pair<std::complex<T>, std::complex<T>> Butterworth<T>::generateBandAnalogPole(int k, T fpw1, T fpw2)
{
    auto thetaK = [pi = PI, order = m_order](int k) -> T {
        return (2 * k - 1) * pi / (2 * order);
    };

    std::complex<T> analogPole(-std::sin(thetaK(k)), std::cos(thetaK(k)));
    std::pair<std::complex<T>, std::complex<T>> poles;
    switch (m_type) {
    case Type::BandReject:
        return poles;
    case Type::BandPass:
    default: {
        std::complex<T> fpw0 = std::sqrt(fpw1 * fpw2);
        std::complex<T> s = T(0.5) * (fpw2 - fpw1) * analogPole / fpw0;
        poles.first = T(2) * PI * fpw0 * (s + std::complex<T>(T(0), T(1)) * std::sqrt(std::complex<T>(T(1), T(0)) - s * s));
        poles.second = T(2) * PI * fpw0 * (s - std::complex<T>(T(0), T(1)) * std::sqrt(std::complex<T>(T(1), T(0)) - s * s));
        return poles;
    }
    }
}

template <typename T>
vectXc_t<T> Butterworth<T>::generateAnalogZeros()
{
    switch (m_type) {
    case Type::HighPass:
        return vectXc_t<T>::Constant(m_order, std::complex<T>(1));
    case Type::BandPass:
        return (vectXc_t<T>(2 * m_order) << vectXc_t<T>::Constant(m_order, std::complex<T>(-1)), vectXc_t<T>::Constant(m_order, std::complex<T>(1))).finished();
    case Type::LowPass:
    default:
        return vectXc_t<T>::Constant(m_order, std::complex<T>(-1));
    }
}

template <typename T>
void Butterworth<T>::scaleAmplitude(Eigen::Ref<vectX_t<T>> aCoeff, Eigen::Ref<vectX_t<T>> bCoeff)
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
    case Type::LowPass:
    default:
        num = aCoeff.sum();
        denum = bCoeff.sum();
        break;
    }

    bCoeff *= num / denum;
}

} // namespace fratio