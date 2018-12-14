#include <limits>

namespace fratio {

// Public static functions
template <typename T>
std::string GenericFilter<T>::filterStatus(FilterStatus status)
{
    switch (status) {
    case FilterStatus::NONE:
        return "Filter is uninitialized";
    case FilterStatus::READY:
        return "Filter is ready to be used";
    case FilterStatus::ALL_COEFF_MISSING:
        return "Filter has none of its coefficient initialized";
    case FilterStatus::A_COEFF_MISSING:
        return "Filter has its 'a' coefficients uninitialized";
    case FilterStatus::A_COEFF_MISSING:
        return "Filter has its 'b' coefficients uninitialized";
    case FilterStatus::BAD_CUTOFF_FREQUENCY:
        return "Filter has a received a bad cut-off frequency. It must be inferior to the sampling frequency";
    default:
        return "I forgot to implement this error documentation";
    }
}

// Public functions

template <typename T>
T GenericFilter<T>::stepFilter(const T& data)
{
    assert(m_status == FilterStatus::READY);

    // Slide data (can't use SIMD, but should be small)
    for (auto rit1 = m_rawData.rbegin(), rit2 = m_rawData.rbegin() + 1; rit2 != m_rawData.rend(); ++rit1, ++rit2)
        *rit1 = *rit2;
    for (auto rit1 = m_filteredData.rbegin(), rit2 = m_filteredData.rbegin() + 1; rit2 != m_filteredData.rend(); ++rit1, ++rit2)
        *rit1 = *rit2;

    m_rawData[0] = data;
    m_filteredData[0] = m_bCoeff.dot(m_rawData) - m_aCoeff.dot(m_filteredData);
    return m_filteredData[0];
}

template <typename T>
Eigen::VectorX<T> GenericFilter<T>::filter(const Eigen::VectorX<T>& data)
{
    Eigen::VectorX<T> results(data.size());
    if (!getFilterResults(results, data))
        return Eigen::VectorX<T>();

    return results;
}

template <typename T>
bool GenericFilter<T>::getFilterResults(Eigen::Ref<Eigen::VectorX<T>> results, const Eigen::VectorX<T>& data)
{
    assert(m_status == FilterStatus::READY);
    if (results.size() != data.size())
        return false;

    T* res = results.data();
    for (T d : data)
        *(res++) = stepFilter(d);

    return true;
}

template <typename T>
void GenericFilter<T>::resetFilter()
{
    m_filteredData.setZero(m_aCoeff.size());
    m_rawData.setZero(m_bCoeff.size());
}

template <typename T>
bool GenericFilter<T>::setCoeffs(const std::vector<T>& aCoeff, const std::vector<T>& bCoeff)
{
    if (!checkCoeffs(aCoeff, bCoeff))
        return false;

    m_aCoeff = Eigen::Map<Eigen::VectorX<T>>(aCoeff.data(), aCoeff.size());
    m_bCoeff = Eigen::Map<Eigen::VectorX<T>>(bCoeff.data(), bCoeff.size());
    resetFilter();
    normalizeCoeffs();
    return true;
}

template <typename T>
void GenericFilter<T>::setCoeffs(const Eigen::VectorX<T>& aCoeff, const Eigen::VectorX<T>& bCoeff)
{
    if (!checkCoeffs(aCoeff, bCoeff))
        return false;

    m_aCoeff = aCoeff;
    m_bCoeff = bCoeff;
    resetFilter();
    normalizeCoeffs();
    return true;
}

template <typename T>
void GenericFilter<T>::getCoeffs(std::vector<T>& aCoeff, std::vector<T>& bCoeff) const
{
    aCoeff.assign(m_aCoeff.data(), m_aCoeff.data() + m_aCoeff.size());
    bCoeff.assign(m_bCoeff.data(), m_bCoeff.data() + m_bCoeff.size());
}

template <typename T>
void GenericFilter<T>::getCoeffs(Eigen::Ref<Eigen::VectorX<T>> aCoeff, Eigen::Ref<Eigen::VectorX<T>> bCoeff) const
{
    aCoeff = m_aCoeff;
    bCoeff = m_bCoeff;
}

// Protected functions

template <typename T>
GenericFilter<T>::GenericFilter(const Eigen::VectorX<T>& aCoeff, const Eigen::VectorX<T>& bCoeff)
    : m_aCoeff(aCoeff)
    , m_bCoeff(bCoeff)
    , m_filteredData(aCoeff.size())
    , m_rawData(bCoeff.size())
{
    if(!checkCoeffs(aCoeff, bCoeff))
        return;

    resetFilter();
    normalizeCoeffs();
}

template <typename T>
void GenericFilter<T>::normalizeCoeffs()
{
    assert(m_status == FilterStatus::READY);

    T a0 = m_aCoeff(0);
    if (std::abs(a0 - T(1)) < std::numeric_limits<T>::epsilon())
        return;

    m_aCoeff /= a0;
    m_bCoeff /= a0;
}

template <typename T>
template <typename T2>
bool GenericFilter<T>::checkCoeffs(const T2& aCoeff, const T2& bCoeff)
{
    using namespace FilterStatus;

    m_status = NONE;
    if (aCoeff.size() == 0)
        m_status = A_COEFF_MISSING;
    else if (std::abs(aCoeff[0]) < std::numeric_limits<T>::epsilon())
        m_status = BAD_A_COEFF;

    if (bCoeff.size() == 0)
        m_status = (m_status == A_COEFF_MISSING ? ALL_COEFF_MISSING : B_COEFF_MISSING);

    if (m_status == NONE)
        m_status = READY;

    return m_status == READY;
}

} // namespace fratio