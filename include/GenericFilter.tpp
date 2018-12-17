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
    case FilterStatus::BAD_ORDER_SIZE:
        return "You try to initialize the filter with an order inferior or equal to 0 (window size for the moving average)";
    case FilterStatus::ALL_COEFF_MISSING:
        return "Filter has none of its coefficient initialized";
    case FilterStatus::A_COEFF_MISSING:
        return "Filter has its 'a' coefficients uninitialized";
    case FilterStatus::B_COEFF_MISSING:
        return "Filter has its 'b' coefficients uninitialized";
    case FilterStatus::BAD_FREQUENCY_VALUE:
        return "Filter has a received a frequency that is negative or equal to zero";
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
    for (Eigen::Index i = m_rawData.size() - 1; i > 0; --i)
        m_rawData(i) = m_rawData(i - 1);
    for (Eigen::Index i = m_filteredData.size() - 1; i > 0; --i)
        m_filteredData(i) = m_filteredData(i - 1);

    m_rawData[0] = data;
    m_filteredData[0] = 0;
    m_filteredData[0] = m_bCoeff.dot(m_rawData) - m_aCoeff.dot(m_filteredData);
    return m_filteredData[0];
}

template <typename T>
vectX_t<T> GenericFilter<T>::filter(const vectX_t<T>& data)
{
    vectX_t<T> results(data.size());
    if (!getFilterResults(results, data))
        return vectX_t<T>();

    return results;
}

template <typename T>
bool GenericFilter<T>::getFilterResults(Eigen::Ref<vectX_t<T>> results, const vectX_t<T>& data)
{
    assert(m_status == FilterStatus::READY);
    if (results.size() != data.size())
        return false;

    for (Eigen::Index i = 0; i < data.size(); ++i)
        results(i) = stepFilter(data(i));

    return true;
}

template <typename T>
void GenericFilter<T>::resetFilter()
{
    m_filteredData.setZero(m_aCoeff.size());
    m_rawData.setZero(m_bCoeff.size());
}

template <typename T>
template <typename T2>
bool GenericFilter<T>::setCoeffs(T2&& aCoeff, T2&& bCoeff)
{
    static_assert(std::is_same_v<T2, vectX_t<T>>, "The coefficents should be of type vectX_t<T>");

    if (!checkCoeffs(aCoeff, bCoeff))
        return false;

    m_aCoeff = aCoeff;
    m_bCoeff = bCoeff;
    resetFilter();
    normalizeCoeffs();
    return true;
}

template <typename T>
void GenericFilter<T>::getCoeffs(vectX_t<T>& aCoeff, vectX_t<T>& bCoeff) const
{
    aCoeff = m_aCoeff;
    bCoeff = m_bCoeff;
}

// Protected functions

template <typename T>
GenericFilter<T>::GenericFilter(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff)
    : m_aCoeff(aCoeff)
    , m_bCoeff(bCoeff)
    , m_filteredData(aCoeff.size())
    , m_rawData(bCoeff.size())
{
    if (!checkCoeffs(aCoeff, bCoeff))
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
bool GenericFilter<T>::checkCoeffs(const vectX_t<T>& aCoeff, const vectX_t<T>& bCoeff)
{
    m_status = FilterStatus::NONE;
    if (aCoeff.size() == 0)
        m_status = FilterStatus::A_COEFF_MISSING;
    else if (std::abs(aCoeff[0]) < std::numeric_limits<T>::epsilon())
        m_status = FilterStatus::BAD_A_COEFF;

    if (bCoeff.size() == 0)
        m_status = (m_status == FilterStatus::A_COEFF_MISSING ? FilterStatus::ALL_COEFF_MISSING : FilterStatus::B_COEFF_MISSING);

    if (m_status == FilterStatus::NONE)
        m_status = FilterStatus::READY;

    return m_status == FilterStatus::READY;
}

} // namespace fratio