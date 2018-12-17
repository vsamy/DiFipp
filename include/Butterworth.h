#pragma once

#include "DigitalFilter.h"
#include "typedefs.h"
#include <cmath>
#include <complex>

namespace fratio {

// https://www.dsprelated.com/showarticle/1119.php
template <typename T>
class Butterworth : public DigitalFilter<T> {
public:
    T PI = static_cast<T>(M_PI);

public:
    enum class Type {
        LowPass,
        HighPass
    };

    // public:
    //     static double minimumRequiredFrequency(...);
public:
    Butterworth(Type type = Type::LowPass);
    Butterworth(int order, T fc, T fs, Type type = Type::LowPass);

    void setFilterParameters(int order, T fc, T fs);

private:
    void initialize(int order, T fc, T fs);
    void computeDigitalRep();
    std::complex<T> generateAnalogPole(T fpw, int k);
    vectX_t<std::complex<T>> generateAnalogZeros();
    void scaleAmplitude(Eigen::Ref<vectX_t<T>> aCoeff, Eigen::Ref<vectX_t<T>> bCoeff);

private:
    Type m_type;
    int m_order;
    T m_fc;
    T m_fs;
    vectX_t<std::complex<T>> m_poles;
};

} // namespace fratio

#include "Butterworth.tpp"
