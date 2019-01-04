#pragma once

#include "type_checks.h"
#include "typedefs.h"

namespace fratio {

/*! \brief Transform an analog signal to a discrete signal and vice versa.
 * 
 * \see https://en.wikipedia.org/wiki/Bilinear_transform
 * \tparam T Floating (complex) types.
 */
template <typename T>
struct BilinearTransform {
    using SubType = internal::complex_sub_type_t<T>; /*!< Sub-type of the complex if T is complex, T otherwise */
    static_assert(std::is_floating_point<SubType>::value, "This struct can only accept floating point types (real and complex).");

    /*! \brief Transformation from analog to discrete.
     * \param fs Sampling frequency.
     * \param sPlanePole Analog data.
     * \param[out] zPlanePole Resulting discrete data.
     */
    static void SToZ(SubType fs, const T& sPlanePole, T& zPlanePole);
    /*! \brief Transformation from analog to discrete.
     * \param fs Sampling frequency.
     * \param sPlanePole Analog signal.
     * \param[out] zPlanePole Resulting discrete signal.
     */
    static void SToZ(SubType fs, const vectX_t<T>& sPlanePoles, Eigen::Ref<vectX_t<T>>& zPlanePoles); // Can be optimized maybe
    /*! \brief Transformation from discrete to analog.
     * \param fs Sampling frequency.
     * \param zPlanePole Discrete data.
     * \param[out] sPlanePole Resulting analog data.
     */
    static void ZToS(SubType fs, const T& zPlanePole, T& sPlanePole);
    /*! \brief Transformation from discrete to analog.
     * \param fs Sampling frequency.
     * \param zPlanePole Discrete signal.
     * \param[out] sPlanePole Resulting analog signal.
     */
    static void ZToS(SubType fs, const vectX_t<T>& zPlanePoles, Eigen::Ref<vectX_t<T>>& sPlanePoles); // Can be optimized maybe
};

template <typename T>
void BilinearTransform<T>::SToZ(SubType fs, const T& sPlanePole, T& zPlanePole)
{
    T scalePole = sPlanePole / (2 * fs);
    zPlanePole = (T(1) + scalePole) / (T(1) - scalePole);
}

template <typename T>
void BilinearTransform<T>::SToZ(SubType fs, const vectX_t<T>& sPlanePoles, Eigen::Ref<vectX_t<T>>& zPlanePoles)
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
void BilinearTransform<T>::ZToS(SubType fs, const vectX_t<T>& zPlanePoles, Eigen::Ref<vectX_t<T>>& sPlanePoles)
{
    assert(zPlanePoles.size() == sPlanePoles.size());
    for (Eigen::Index k = 0; k < sPlanePoles.size(); ++k)
        ZToS(fs, zPlanePoles(k), sPlanePoles(k));
}

} // namespace fratio