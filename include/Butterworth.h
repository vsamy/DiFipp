#pragma once

#include "GenericFilter.h"
#include <complex>
#include <vector>

namespace fratio {

// https://www.dsprelated.com/showarticle/1119.php
template <typename T>
class Butterworth : public GenericFilter<T> {
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
    void transformFilter();

private:
    Type m_type;
    size_t m_order;
    T m_fc;
    T m_fs;
    std::vector<std::complex<T>> m_poles;
};

} // namespace fratio

#include "Butterworth.inl"