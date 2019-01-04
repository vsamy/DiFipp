#pragma once

#include "DigitalFilter.h"
#include "typedefs.h"
#include <cmath>
#include <complex>

namespace fratio {

// https://www.dsprelated.com/showarticle/1119.php
// https://www.dsprelated.com/showarticle/1135.php
// https://www.dsprelated.com/showarticle/1128.php
// https://www.dsprelated.com/showarticle/1131.php
// https://www.mathworks.com/help/signal/ref/butter.html
template <typename T>
class Butterworth : public DigitalFilter<T> {
public:
    static T PI;

public:
    enum class Type {
        LowPass,
        HighPass,
        BandPass,
        BandReject
    };

public:
    // http://www.matheonics.com/Tutorials/Butterworth.html#Paragraph_3.2
    // https://www.mathworks.com/help/signal/ref/buttord.html#d120e11079
    static std::pair<int, T> findMinimumButter(T wPass, T wStop, T APass, T AStop); // Only for low and high pass

public:
    Butterworth(Type type = Type::LowPass);
    Butterworth(int order, T fc, T fs, Type type = Type::LowPass);
    Butterworth(int order, T fLower, T fUpper, T fs, Type type = Type::BandPass);

    void setFilterParameters(int order, T fc, T fs);
    void setFilterParameters(int order, T fLower, T fUpper, T fs);

private:
    void initialize(int order, T f1, T f2, T fs); // f1 = fc for LowPass/HighPass filter
    void computeDigitalRep(T fc);
    void computeBandDigitalRep(T fLower, T fUpper);
    std::complex<T> generateAnalogPole(int k, T fpw);
    std::pair<std::complex<T>, std::complex<T>> generateBandAnalogPole(int k, T fpw0, T bw);
    vectXc_t<T> generateAnalogZeros(T f0 = T());
    void scaleAmplitude(const vectX_t<T>& aCoeff, Eigen::Ref<vectX_t<T>> bCoeff, const std::complex<T>& bpS = T());

private:
    Type m_type;
    int m_order;
    T m_fs;
    vectXc_t<T> m_poles;
};

} // namespace fratio

#include "Butterworth.tpp"
