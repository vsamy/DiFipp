#include "msfc/filter/GenericFilter.h"
#include "msfc/logging.h"

namespace msfc {

namespace filt {

    GenericFilter::GenericFilter(size_t nData)
    // : m_filteredData(nData, order)
    // , m_rawData(nData, order)
    {
    }

    GenericFilter::GenericFilter(size_t nData, const std::vector<double>& aCoeff, const std::vector<double>& bCoeff)
        : m_aCoeff(aCoeff)
        , m_bCoeff(bCoeff)
        , m_filteredData(aCoeff.size())
        , m_rawData(bCoeff.size())
    {
        assert(aCoeff.size() > 0);
        assert(bCoeff.size() > 0);
        normalize();
    }

    void GenericFilter::setNData(size_t nData)
    {
        // m_filteredData.resize(nData, m_filteredData.cols());
        // m_rawData.resize(nData, nData.cols());
    }

    void GenericFilter::setCoeff(const std::vector<double>& aCoeff, const std::vector<double>& bCoeff)
    {
        assert(aCoeff.size() > 0);
        assert(bCoeff.size() > 0);

        m_nACoeffFilteredData = 0;
        m_nBCoeffFilteredData = 0;
        m_aCoeff = aCoeff;
        m_bCoeff = bCoeff;
        m_filteredData.resize(aCoeff.size());
        m_rawData.resize(bCoeff.size());

        normalize();
    }

    void GenericFilter::getCoeff(std::vector<double>& aCoeff, std::vector<double>& bCoeff) const noexcept
    {
        aCoeff = m_aCoeff;
        bCoeff = m_bCoeff;
    }

    // bool GenericFilter::stepFilter(const Eigen::VectorXd& data)
    // {
    //     if (m_filteredData.rows() != data.size()) {
    //         LOG_ERROR("Bad data size. Expected vector of size " << m_filteredData.rows() << ", got a vector of size " << data.size());
    //         return false;
    //     }

    //     if (m_nFilteredData == 0) {
    //         m_rawData.row(0) = data;
    //         m_filteredData.row(0).noalias() = m_bCoeff(0) * data;
    //         ++m_nFilteredData;
    //     } else if (m_nFilteredData < m_order) {
    //         m_filteredData.row(m_nFilteredData).noalias() = m_bCoeff(0) * data;
    //         for (size_t i = 1; i <= m_nFilteredData; ++i) {
    //             m_filteredData.row(m_nFilteredData).noalias() += m_bCoeff(i) * m_rawData.row(i - 1);
    //         }
    //         ++m_nFilteredData;
    //     }

    //     if (m_nFilteredData >= m_order) {

    //     } else {
    //         for (size_t i = 0; i < m_nFilteredData; ++i) {
    //         }

    //         ++m_nFilteredData;
    //     }
    // }

    // https://stackoverflow.com/questions/50511549/meaning-of-rational-transfer-function-underlying-matlab-filter-or-scipy-signal-f
    bool GenericFilter::stepFilter(double data)
    {
        double filtData;
        for (size_t i = 0; i < m_nBCoeffFilteredData; ++i)
            filtData += m_bCoeff[i] * m_rawData[m_nBCoeffFilteredData - i];
        for (size_t i = 1; i < m_nACoeffFilteredData; ++i)
            filtData -= m_aCoeff[i] * m_filteredData[m_nACoeffFilteredData - i];

        m_filteredData[m_nACoeffFilteredData] = filtData;
        m_rawData[m_nBCoeffFilteredData] = data;
        ++m_nACoeffFilteredData;
        ++m_nBCoeffFilteredData;

        if (m_nACoeffFilteredData == m_filteredData.size()) {
            double* fd = m_filteredData.data();
            for (size_t i = 0; i < m_filteredData.size() - 1; ++i)
                *(fd) = *(++fd);

            --m_nACoeffFilteredData;
        }

        if (m_nBCoeffFilteredData == m_rawData.size()) {
            double* rd = m_rawData.data();
            for (size_t i = 0; i < m_rawData.size() - 1; ++i)
                *(rd) = *(++rd);

            --m_nBCoeffFilteredData;
        }
    }

    void GenericFilter::normalize()
    {
        double a0 = m_aCoeff.front();
        if (std::abs(a0) < 1e-6) // Divide by zero
            LOG_ERROR_AND_THROW(std::invalid_argument, "By filtering value for coefficient a0. Should be superior to 1e-6");

        for (double& a : m_aCoeff)
            a /= a0;
        for (double& b : m_bCoeff)
            b /= a0;
    }

} // namespace filt

} // namespace msfc