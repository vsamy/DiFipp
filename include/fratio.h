#pragma once

#include "Butterworth.h"
#include "DigitalFilter.h"
#include "MovingAverage.h"
#include "polynome_functions.h"

namespace fratio {

using DigitalFilterf = DigitalFilter<float>;
using DigitalFilterd = DigitalFilter<double>;
using MovingAveragef = MovingAverage<float>;
using MovingAveraged = MovingAverage<double>;
using Butterworthf = Butterworth<float>;
using Butterworthd = Butterworth<double>;

using VietaAlgof = VietaAlgo<float>;
using VietaAlgod = VietaAlgo<double>;
using VietaAlgoi = VietaAlgo<int>;
using VietaAlgocf = VietaAlgo<std::complex<float>>;
using VietaAlgocd = VietaAlgo<std::complex<double>>;
using VietaAlgoci = VietaAlgo<std::complex<int>>;

} // namespace fratio