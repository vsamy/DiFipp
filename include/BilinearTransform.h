#pragma once

#include "type_checks.h"
#include "typedefs.h"

namespace fratio {

template <typename T>
struct BilinearTransform {
    using SubType = internal::complex_sub_type_t<T>;
    static_assert(std::is_floating_point<SubType>::value, "This struct can only accept floating point types (real and complex).");

    static void SToZ(SubType fs, const T& sPlanePole, T& zPlanePole);
    static void SToZ(SubType fs, const Eigen::VectorX<T>& sPlanePoles, Eigen::Ref<Eigen::VectorX<T>>& zPlanePoles); // Can be optimized
    static void ZToS(SubType fs, const T& zPlanePole, T& sPlanePole);
    static void ZToS(SubType fs, const Eigen::VectorX<T>& zPlanePoles, Eigen::Ref<Eigen::VectorX<T>>& sPlanePoles); // Can be optimized
};

template <typename T>
void BilinearTransform<T>::SToZ(SubType fs, const T& sPlanePole, T& zPlanePole)
{
    T scalePole = sPlanePole / (2 * fs);
    zPlanePole = (T(1) + scalePole) / (T(1) - scalePole);
}

template <typename T>
void BilinearTransform<T>::SToZ(SubType fs, const Eigen::VectorX<T>& sPlanePoles, Eigen::Ref<Eigen::VectorX<T>>& zPlanePoles)
{
    assert(sPlanePoles.size() == zPlanePoles.size());
    for (Eigen::Index k = 0; k < sPlanePoles.size(); ++k)
        SToZ(fs, sPlanePoles(k), zPlanePoles(k));
}

template <typename T>
void BilinearTransform<T>::ZToS(SubType fs, const T& zPlanePole, T& sPlanePole)
{
    T invPole = T(1) / zPlanePole;
    sPlanePole = 2 * fs * (T(1) - invPole) / (T(1) + invPole);
}

template <typename T>
void BilinearTransform<T>::ZToS(SubType fs, const Eigen::VectorX<T>& zPlanePoles, Eigen::Ref<Eigen::VectorX<T>>& sPlanePoles)
{
    assert(zPlanePoles.size() == sPlanePoles.size());
    for (Eigen::Index k = 0; k < sPlanePoles.size(); ++k)
        ZToS(fs, zPlanePoles(k), sPlanePoles(k));
}

} // namespace fratio