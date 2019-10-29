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

#pragma once

#include "DigitalFilter.h"
#include "typedefs.h"
#include <complex>

namespace difi {

/*! \brief Butterworth digital filter.
 * 
 * This is a specialization of a digital filter in order to use a Butterworth filter.
 * \see https://www.dsprelated.com/showarticle/1119.php
 * \see https://www.dsprelated.com/showarticle/1135.php
 * \see https://www.dsprelated.com/showarticle/1128.php
 * \see https://www.dsprelated.com/showarticle/1131.php
 * \see https://www.mathworks.com/help/signal/ref/butter.html
 * \tparam Floating type.
 */
template <typename T>
class Butterworth : public DigitalFilter<T> {
public:
    /*! \brief Type of butterworth filter0 */
    enum class Type {
        LowPass, /*!< Define a low-pass filter */
        HighPass, /*!< Define a high-pass filter */
        BandPass, /*!< Define a band-pass filter */
        BandReject /*!< Define a band-reject filter */
        // AllPass
        // CombPass
    };

public:
    /*! \brief Function to help you design a Butterworth filter.
     * 
     * It finds optimal values of the order and cut-off frequency.
     * \note Works only for low-pass and high-pass filters.
     * \see http://www.matheonics.com/Tutorials/Butterworth.html#Paragraph_3.2
     * \see https://www.mathworks.com/help/signal/ref/buttord.html#d120e11079
     * \param wPass Pass band edge.
     * \param wStop Stop band edge.
     * \param APass Maximum pass band attenuation.
     * \param AStop Minimum stop band attenuation.
     * \return A pair of { filter order, cut-off frequency }.
     */
    static std::pair<int, T> findMinimumButter(T wPass, T wStop, T APass, T AStop); // Only for low and high pass

public:
    /*! \brief Uninitialized constructor. 
     * \param type Filter type. Default is LowPass.
     */
    Butterworth(Type type = Type::LowPass);
    /*! \brief Constructor for both low-pass and high-pass filters.
     * \param order Order of the filter.
     * \param fc Cut-off frequency.
     * \param fs Sampling frequency.
     * \param type Filter type. Default is LowPass.
     */
    Butterworth(int order, T fc, T fs, Type type = Type::LowPass);
    /*! \brief Constructor for both band-pass and band-reject filters.
     * \param order Order of the filter.
     * \param fLower Lower bound frequency.
     * \param fUpper Upper bound frequency.
     * \param fs Sampling frequency.
     * \param type Filter type. Default is BandPass.
     */
    Butterworth(int order, T fLower, T fUpper, T fs, Type type = Type::BandPass);
    /*! \brief Set filter set of parameters.
     * \param order Order of the filter.
     * \param fc Cut-off frequency.
     * \param fs Sampling frequency.
     */
    void setFilterParameters(int order, T fc, T fs);
    /*! \brief Set filter set of parameters.
     * \param order Order of the filter.
     * \param fLower Lower bound frequency.
     * \param fUpper Upper bound frequency.
     * \param fs Sampling frequency.
     */
    void setFilterParameters(int order, T fLower, T fUpper, T fs);

private:
    /*! \brief Initialize the filter.
     * \param order Order of the filter.
     * \param f1 First frequency parameter.
     * \param f2 Second frequency parameter.
     * \param fs Sampling frequency.
     */
    void initialize(int order, T f1, T f2, T fs);
    /*! \brief Compute the digital filter representation for low-pass and high-pass.
     * \param fc Cut-off frequency.
     */
    void computeDigitalRep(T fc);
    /*! \brief Compute the digital filter representation for band-pass and band-reject.
     * \param fLower Lower bound frequency.
     * \param fUpper Upper bound frequency.
     */
    void computeBandDigitalRep(T fLower, T fUpper);
    /*! \brief Generate an analog pole on the unit circle for low-pass and high-pass.
     * \param k Step on the unit circle.
     * \param fpw Continuous pre-warp cut-off frequency.
     * \return Generated pole.
     */
    std::complex<T> generateAnalogPole(int k, T fpw);
    /*! \brief Generate an analog pole on the unit circle for band-pass and band-reject.
     * \param k Step on the unit circle.
     * \param fpw0 Continuous pre-warp frequency at geometric center.
     * \param bw Bandwith.
     * \return Pair of generated pole.
     */
    std::pair<std::complex<T>, std::complex<T>> generateBandAnalogPole(int k, T fpw0, T bw);
    /*! \brief Generate all analog zeros.
     * \param fpw0 Continuous pre-warp frequency at geometric center (Only use by the band-reject).
     * \return Set of generated zeros.
     */
    vectXc_t<T> generateAnalogZeros(T fpw0 = T());
    /*! \brief Scale coefficients.
     * \param aCoeff Unscaled poles.
     * \param bCoeff Unscaled zeros. 
     * \param bpS Result of \f$H(fc)=1\f$ with \f$H(z)\f$ the discrete transfer function at the geometric center \f$fc\f$ (Only use by band-pass).
     */
    void scaleAmplitude(const vectX_t<T>& aCoeff, Eigen::Ref<vectX_t<T>> bCoeff, const std::complex<T>& bpS = T());

private:
    Type m_type; /*!< Filter type */
    int m_order; /*!< Filter order */
    T m_fs; /*!< Filter sampling frequency */
};

} // namespace difi

#include "Butterworth.tpp"
