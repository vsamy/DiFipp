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

#include "BilinearTransform.h"
#include "Butterworth.h"
#include "DigitalFilter.h"
#include "GenericFilter.h"
#include "MovingAverage.h"
#include "differentiators.h"
#include "polynome_functions.h"
#include "typedefs.h"

namespace difi {

// Filters
using DigitalFilterf = DigitalFilter<float>;
using DigitalFilterd = DigitalFilter<double>;
using MovingAveragef = MovingAverage<float>;
using MovingAveraged = MovingAverage<double>;
using Butterworthf = Butterworth<float>;
using Butterworthd = Butterworth<double>;

// Polynome helper functions
using VietaAlgof = VietaAlgo<float>;
using VietaAlgod = VietaAlgo<double>;
using VietaAlgoi = VietaAlgo<int>;
using VietaAlgocf = VietaAlgo<std::complex<float>>;
using VietaAlgocd = VietaAlgo<std::complex<double>>;
using VietaAlgoci = VietaAlgo<std::complex<int>>;

// Bilinear transformation functions
using BilinearTransformf = BilinearTransform<float>;
using BilinearTransformd = BilinearTransform<double>;
using BilinearTransformcf = BilinearTransform<std::complex<float>>;
using BilinearTransformcd = BilinearTransform<std::complex<double>>;

// 1st order centered differentiators
template <int N> using CenteredDiffBasicf = CenteredDiffBasic<float, N>;
template <int N> using CenteredDiffBasicd = CenteredDiffBasic<double, N>;
template <int N> using CenteredDiffLowNoiseLanczosf = CenteredDiffLowNoiseLanczos<float, N>;
template <int N> using CenteredDiffLowNoiseLanczosd = CenteredDiffLowNoiseLanczos<double, N>;
template <int N> using CenteredDiffSuperLowNoiseLanczosf = CenteredDiffSuperLowNoiseLanczos<float, N>;
template <int N> using CenteredDiffSuperLowNoiseLanczosd = CenteredDiffSuperLowNoiseLanczos<double, N>;
template <int N> using CenteredDiffNoiseRobust2f = CenteredDiffNoiseRobust2<float, N>;
template <int N> using CenteredDiffNoiseRobust2d = CenteredDiffNoiseRobust2<double, N>;
template <int N> using CenteredDiffNoiseRobust4f = CenteredDiffNoiseRobust4<float, N>;
template <int N> using CenteredDiffNoiseRobust4d = CenteredDiffNoiseRobust4<double, N>;
template <int N> using TVCenteredDiffNoiseRobust2f = TVCenteredDiffNoiseRobust2<float, N>;
template <int N> using TVCenteredDiffNoiseRobust2d = TVCenteredDiffNoiseRobust2<double, N>;
template <int N> using TVCenteredDiffNoiseRobust4f = TVCenteredDiffNoiseRobust4<float, N>;
template <int N> using TVCenteredDiffNoiseRobust4d = TVCenteredDiffNoiseRobust4<double, N>;

// 2nd order centered differentiators
template <int N> using CenteredDiffSecondOrderf = CenteredDiffSecondOrder<float, N>;
template <int N> using CenteredDiffSecondOrderd = CenteredDiffSecondOrder<double, N>;
template <int N> using TVCenteredDiffSecondOrderf = TVCenteredDiffSecondOrder<float, N>;
template <int N> using TVCenteredDiffSecondOrderd = TVCenteredDiffSecondOrder<double, N>;

// 1st order backward differentiators
template <int N> using BackwardDiffNoiseRobustf = BackwardDiffNoiseRobust<float, N>;
template <int N> using BackwardDiffNoiseRobustd = BackwardDiffNoiseRobust<double, N>;
template <int N> using BackwardDiffHybridNoiseRobustf = BackwardDiffHybridNoiseRobust<float, N>;
template <int N> using BackwardDiffHybridNoiseRobustd = BackwardDiffHybridNoiseRobust<double, N>;
template <int N> using TVBackwardDiffNoiseRobustf = TVBackwardDiffNoiseRobust<float, N>;
template <int N> using TVBackwardDiffNoiseRobustd = TVBackwardDiffNoiseRobust<double, N>;
template <int N> using TVBackwardDiffHybridNoiseRobustf = TVBackwardDiffHybridNoiseRobust<float, N>;
template <int N> using TVBackwardDiffHybridNoiseRobustd = TVBackwardDiffHybridNoiseRobust<double, N>;

// 2nd order backward differentiators
template <int N> using BackwardDiffSecondOrderf = BackwardDiffSecondOrder<float, N>;
template <int N> using BackwardDiffSecondOrderd = BackwardDiffSecondOrder<double, N>;
template <int N> using TVBackwardDiffSecondOrderf = TVBackwardDiffSecondOrder<float, N>;
template <int N> using TVBackwardDiffSecondOrderd = TVBackwardDiffSecondOrder<double, N>;

} // namespace difi