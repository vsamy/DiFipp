#pragma once

#include <stddef.h>
#include <vector>

namespace fratio {

class GenericFilter {
public:
    GenericFilter() = default;
    GenericFilter(size_t nData);
    GenericFilter(size_t nData, const std::vector<double>& aCoeff, const std::vector<double>& bCoeff);

    void setNData(size_t nData);
    void setCoeff(const std::vector<double>& aCoeff, const std::vector<double>& bCoeff);
    void getCoeff(std::vector<double>& aCoeff, std::vector<double>& bCoeff) const noexcept;
    size_t aOrder() const noexcept { return m_aCoeff.size(); }
    size_t bOrder() const noexcept { return m_bCoeff.size(); }

    // bool stepFilter(const Eigen::VectorXd& data);
    // https://stackoverflow.com/questions/50511549/meaning-of-rational-transfer-function-underlying-matlab-filter-or-scipy-signal-f
    double stepFilter(double data);
    std::vector<double> filter(const std::vector<double>& data);
    void resetFilter();

    // Eigen::VectorXd results() const noexcept;
    // std::vector<double> results() const noexcept { return m_filteredData; }

private:
    void checkCoeffs(const std::vector<double>& aCoeff, const std::vector<double>& bCoeff);
    void normalize();
    void shiftData();

protected:
    std::vector<double> m_aCoeff;
    std::vector<double> m_bCoeff;

    std::vector<double> m_filteredData;
    std::vector<double> m_rawData;
    // Eigen::MatrixXd m_filteredData;
    // Eigen::MatrixXd m_rawData;
};

} // namespace fratio