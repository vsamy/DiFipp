#include "GenericFilter.h"

#include <sstream>

namespace fratio {

GenericFilter::GenericFilter(size_t nData)
{
}

GenericFilter::GenericFilter(size_t nData, const std::vector<double>& aCoeff, const std::vector<double>& bCoeff)
    : m_aCoeff(aCoeff)
    , m_bCoeff(bCoeff)
    , m_filteredData(aCoeff.size(), 0)
    , m_rawData(bCoeff.size(), 0)
{
    checkCoeffs(aCoeff, bCoeff);
    normalize();
}

void GenericFilter::setNData(size_t nData)
{
    // m_filteredData.resize(nData, m_filteredData.cols());
    // m_rawData.resize(nData, nData.cols());
}

void GenericFilter::setCoeff(const std::vector<double>& aCoeff, const std::vector<double>& bCoeff)
{
    checkCoeffs(aCoeff, bCoeff);

    m_aCoeff = aCoeff;
    m_bCoeff = bCoeff;
    resetFilter();

    normalize();
}

void GenericFilter::getCoeff(std::vector<double>& aCoeff, std::vector<double>& bCoeff) const noexcept
{
    aCoeff = m_aCoeff;
    bCoeff = m_bCoeff;
}

double GenericFilter::stepFilter(double data)
{
    // Slide data
    for (auto rit1 = m_rawData.rbegin(), rit2 = m_rawData.rbegin() + 1; rit2 != m_rawData.rend(); ++rit1, ++rit2)
        *rit1 = *rit2;
    for (auto rit1 = m_filteredData.rbegin(), rit2 = m_filteredData.rbegin() + 1; rit2 != m_filteredData.rend(); ++rit1, ++rit2)
        *rit1 = *rit2;

    double filtData = 0.;
    m_rawData[0] = data;

    for (size_t k = 0; k < m_bCoeff.size(); ++k)
        filtData += m_bCoeff[k] * m_rawData[k];
    for (size_t k = 1; k < m_aCoeff.size(); ++k)
        filtData -= m_aCoeff[k] * m_filteredData[k];

    m_filteredData[0] = filtData;

    return filtData;
}

std::vector<double> GenericFilter::filter(const std::vector<double>& data)
{
    std::vector<double> results;
    results.reserve(data.size());

    for (double d : data)
        results.emplace_back(stepFilter(d));

    return results;
}

void GenericFilter::resetFilter()
{
    m_filteredData.assign(m_aCoeff.size(), 0);
    m_rawData.assign(m_bCoeff.size(), 0);
}

void GenericFilter::normalize()
{
    double a0 = m_aCoeff.front();
    if (std::abs(a0) < 1e-6) // Divide by zero
        throw std::invalid_argument("By filtering value for coefficient a0. Should be superior to 1e-6");

    for (double& a : m_aCoeff)
        a /= a0;
    for (double& b : m_bCoeff)
        b /= a0;
}

void GenericFilter::checkCoeffs(const std::vector<double>& aCoeff, const std::vector<double>& bCoeff)
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