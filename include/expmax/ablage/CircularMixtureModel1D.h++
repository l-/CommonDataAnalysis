///**
// * @file CircularMixtureModel1D.h++
// *
// * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
// *
// *  Created on: Jan 5, 2011
// *
// */
//
//#pragma once
//
//#include "math/CircularStatistics.h++"
//#include "FitMultivariateMulticlassByEM.h++"
//
//namespace CDA {
//
//
///**
// * @class CircularMixtureModel1D
// *
// * A univariate VonMises-MixtureModel class, with EM fitting routine built in
// *
// */
//class CircularMixtureModel1D : public FitUnivariateMulticlassByEM {
//
//    typedef double datapoint_t;
//
//    using EM<datapoint_t>::m_theta;
//    using EM<datapoint_t>::m_data;
//
//    using FitUnivariateMulticlassByEM::K;
//
//protected:
//
//    /**
//     * @brief <b>E-step</b>: update hidden attributes, i.e. class membership probabilities
//     */
//    void update_thetas();
//
//public:
//
//    /**
//     * @brief Evaluate the cluster model PDF with class parameters \f$k\f$ at \f$\vec{x}\f$
//     *
//     * @param[in] k class
//     * @param[in] x feature vector
//     *
//     * @return probability \f$p(x,k\vert\theta)\f$
//     */
//    double evalPDF(const unsigned k, const datapoint_t& x) const;
//
//    /**
//     * @brief Call superclass constructor to fill data fields
//     */
//    CircularMixtureModel1D(const unsigned K_)
//    : FitUnivariateMulticlassByEM(K_)
//    {  }
//
//    /**
//     * @brief Estimated parameter getter.
//     * Get current estimated mean vector of class k
//     *
//     * @param[in] k class no.
//     */
//    double getMu(const unsigned k) const;
//
//    /**
//     * @brief Estimated parameter getter
//     * Get current estimated sigma of class k
//     *
//     * @param[in] k class no.
//     */
//    double getKappa(const unsigned k) const;
//
//    /**
//     * Evaluate the PDF with parameters of class k
//     *
//     * @param[in] k class no.
//     * @param[in] x
//     */
//    inline double evalPDF(const unsigned k, const double x) const;
//
//    /**
//     * @brief Initialize known samples, i.e. your data points
//     */
//    template<class II>
//    void setData(std::pair<II, II> data_) {
//        FitUnivariateMulticlassByEM::setData(data_);
//    }
//
//    /**
//     * @brief Class name for logging ...
//     */
//    const std::string className() const;
//
//    /**
//     * @brief Each class gets its own line, that's the easiest way
//     */
//    const std::string getCSVHeader() const;
//
//};
//
//
//} //namespace
