#pragma once

#include "DigitalFilter.h"
#include "typedefs.h"
#include <cmath>
#include <complex>

namespace fratio {

// https://www.dsprelated.com/showarticle/1119.php
// https://www.dsprelated.com/showarticle/1135.php
// https://www.dsprelated.com/showarticle/1128.php
template <typename T>
class Butterworth : public DigitalFilter<T> {
public:
    T PI = static_cast<T>(M_PI);

public:
    enum class Type {
        LowPass,
        HighPass,
        BandPass,
        BandReject
    };

    // public:
    //     static double minimumRequiredFrequency(...);
public:
    Butterworth(Type type = Type::LowPass);
    Butterworth(int order, T fc, T fs, Type type = Type::LowPass);
    Butterworth(int order, T bw, T fs, T fCenter, Type type = Type::BandPass);

    void setFilterParameters(int order, T fc, T fs, T fCenter = T(0));

private:
    void initialize(int order, T fc, T fs, T fCenter = T(0)); // fc = bw for bandPass filter
    void computeDigitalRep(T fc);
    void computeBandDigitalRep(T bw, T fCenter);
    std::complex<T> generateAnalogPole(int k, T fpw);
    std::pair<std::complex<T>, std::complex<T>> generateBandAnalogPole(int k, T fpw1, T fpw2);
    vectXc_t<T> generateAnalogZeros();
    void scaleAmplitude(Eigen::Ref<vectX_t<T>> aCoeff, Eigen::Ref<vectX_t<T>> bCoeff);

private:
    Type m_type;
    int m_order;
    T m_fc;
    T m_bw;
    T m_fs;
    T m_fCenter;
    vectXc_t<T> m_poles;
};

} // namespace fratio

#include "Butterworth.tpp"
