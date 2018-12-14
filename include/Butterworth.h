#pragma once

#include "GenericFilter.h"
#include "typedefs.h"
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
    Butterworth(size_t order, T fc, T fs, Type type = Type::LowPass);

    void setFilterParameters(size_t order, T fc, T fs);

private:
    void initialize(size_t order, T fc, T fs);
    void computeDigitalRep();
    void updateCoeffSize();
    std::complex<T> generateAnalogPole(T fpw, size_t k);
    std::vector<std::complex<T>> generateAnalogZeros();
    void scaleAmplitude();

private:
    Type m_type;
    size_t m_order;
    T m_fc;
    T m_fs;
    Eigen::VectorX<std::complex<T>> m_poles;
};

} // namespace fratio

#include "Butterworth.tpp"
