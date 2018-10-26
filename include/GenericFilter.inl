namespace fratio {

// Public functions

template <typename T>
T GenericFilter<T>::stepFilter(T data)
{
    // Slide data
    for (auto rit1 = m_rawData.rbegin(), rit2 = m_rawData.rbegin() + 1; rit2 != m_rawData.rend(); ++rit1, ++rit2)
        *rit1 = *rit2;
    for (auto rit1 = m_filteredData.rbegin(), rit2 = m_filteredData.rbegin() + 1; rit2 != m_filteredData.rend(); ++rit1, ++rit2)
        *rit1 = *rit2;

    T filtData = 0.;
    m_rawData[0] = data;

    for (size_t k = 0; k < m_bCoeff.size(); ++k)
        filtData += m_bCoeff[k] * m_rawData[k];
    for (size_t k = 1; k < m_aCoeff.size(); ++k)
        filtData -= m_aCoeff[k] * m_filteredData[k];

    m_filteredData[0] = filtData;

    return filtData;
}

template <typename T>
std::vector<T> GenericFilter<T>::filter(const std::vector<T>& data)
{
    std::vector<T> results;
    results.reserve(data.size());

    for (T d : data)
        results.emplace_back(stepFilter(d));

    return results;
}

template <typename T>
void GenericFilter<T>::resetFilter()
{
    m_filteredData.assign(m_aCoeff.size(), 0);
    m_rawData.assign(m_bCoeff.size(), 0);
}

template <typename T>
void GenericFilter<T>::getCoeff(std::vector<T>& aCoeff, std::vector<T>& bCoeff) const noexcept
{
    aCoeff = m_aCoeff;
    bCoeff = m_bCoeff;
}

// Protected functions

template <typename T>
GenericFilter<T>::GenericFilter(const std::vector<T>& aCoeff, const std::vector<T>& bCoeff)
    : m_aCoeff(aCoeff)
    , m_bCoeff(bCoeff)
    , m_filteredData(aCoeff.size(), 0)
    , m_rawData(bCoeff.size(), 0)
{
    checkCoeff(aCoeff, bCoeff);
    normalize();
}

template <typename T>
void GenericFilter<T>::normalize()
{
    T a0 = m_aCoeff.front();
    if (std::abs(a0) < 1e-8) // Divide by zero
        throw std::invalid_argument("By filtering value for coefficient a0. Should be superior to 1e-8");

    if (std::abs(a0 - 1) < 1e-8)
        return;

    for (T& a : m_aCoeff)
        a /= a0;
    for (T& b : m_bCoeff)
        b /= a0;
}

template <typename T>
void GenericFilter<T>::checkCoeff(const std::vector<T>& aCoeff, const std::vector<T>& bCoeff)
{
    std::stringstream err;
    if (aCoeff.size() == 0)
        err << "The size of coefficient 'a' should greater than 0\n";
    if (bCoeff.size() == 0)
        err << "The size of coefficient 'b' should greater than 0\n";

    if (err.str().size() > 0)
        throw std::runtime_error(err.str());
}

} // namespace fratio