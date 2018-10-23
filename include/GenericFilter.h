#pragma once

#include <Eigen/Core>
#include <stddef.h>
#include <vector>

namespace msfc {

namespace filt {

    class GenericFilter {
    public:
        GenericFilter() = default;
        GenericFilter(size_t nData);
        GenericFilter(size_t nData, const std::vector<double>& aCoeff, const std::vector<double>& bCoeff);

        void setNData(size_t nData);
        void setCoeff(const std::vector<double>& aCoeff, const std::vector<double>& bCoeff);
        void getCoeff(std::vector<double>& aCoeff, std::vector<double>& bCoeff) const noexcept;
        size_t filterOrder() const noexcept { return m_order; }

        // bool stepFilter(const Eigen::VectorXd& data);
        bool stepFilter(double data);

        // Eigen::VectorXd results() const noexcept;
        double results() const noexcept;

    private:
        void normalize();
        void shiftData();

    protected:
        size_t m_nACoeffFilteredData;
        size_t m_nBCoeffFilteredData;
        std::vector<double> m_aCoeff;
        std::vector<double> m_bCoeff;

        std::vector<double> m_filteredData;
        std::vector<double> m_rawData;
        // Eigen::MatrixXd m_filteredData;
        // Eigen::MatrixXd m_rawData;
    };

} // namespace filt

} // namespace msfc