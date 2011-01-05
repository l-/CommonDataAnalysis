/**
 * @file GaussianMixtureModel1D.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 5, 2011
 *
 * @section NOTE
 * This one, the easiest one, actually works.
 *
 */



#pragma once

#include "FitMultivariateMulticlassByEM.h++"


namespace CDA {


/**
 * @class GaussianMixtureModel1D
 *
 * A GaussianMixtureModel class, with EM fitting routine built in
 *
 */
class GaussianMixtureModel1D : public FitUnivariateMulticlassByEM {

public:
    typedef double datapoint_t;

private:
    using EM<datapoint_t>::m_theta;
    using EM<datapoint_t>::m_data;

    using FitUnivariateMulticlassByEM::K;

    inline double squaredDistanceToMean(const unsigned k, const datapoint_t x) const {
        return fabs(x - getMean(k)); // OBS! fabs!
    }

protected:

    /**
     * @brief <b>E-step</b>: update hidden attributes, i.e. class membership probabilities
     */
    void update_thetas();

public:

    /**
     * @brief Evaluate the cluster model PDF with class parameters \f$k\f$ at \f$\vec{x}\f$
     *
     * @param[in] k class
     * @param[in] x feature vector
     *
     * @return probability \f$p(x,k\vert\theta)\f$
     */
    inline double evalPDF(const unsigned k, const datapoint_t& x) const;

    /**
     * @brief Call superclass constructor to fill data fields
     */
    GaussianMixtureModel1D(const unsigned K_)
    : FitUnivariateMulticlassByEM(K_)
    {  }

    /**
     * @brief Estimated parameter getter.
     * Get current estimated mean vector of class k
     *
     * @param[in] k class no.
     */
    double getMean(const unsigned k) const;

    /**
     * @brief Estimated parameter getter
     * Get current estimated sigma of class k
     *
     * @param[in] k class no.
     */
    double getSigma(const unsigned k) const;
    //
    //    /**
    //     * Evaluate the PDF with parameters of class k
    //     *
    //     * @param[in] k class no.
    //     * @param[in] x
    //     */
    //    inline double evalPDF(const unsigned k, const double x) const;

    /**
     * @brief Initialize known samples, i.e. your data points
     *
     * @param[in] data_ A pair of iterators.
     * @todo we should check II for compliance.
     */
    template<class II>
    void setData(std::pair<II, II> data_) {
        FitUnivariateMulticlassByEM::setData(data_);
    }

    /**
     * @brief Find discriminant between two classes
     *
     * @param[in] k1 a class
     * @param[in] k2 a class
     *
     * At least this. Is there an easier way?
     */
    boost::function<double(const double)>
    getDiscriminantFunction(const int k1, const int k2) const {
        // @todo double-check

        assert(k1 != k2);

        return boost::function<double(const double)>(
                boost::lambda::bind(fabs, (boost::lambda::_1 - getMean(k1))) * pow(getSigma(k2),2)
                - boost::lambda::bind(fabs, (boost::lambda::_1 - getMean(k2)))  * pow(getSigma(k1),2)
                - 2 * log ( getSigma(k2) / getSigma(k1) ) * pow(getSigma(k2),2) * pow(getSigma(k1),2));
        //          < 0  ? k1 : k2 );

        // eh? ist auch falsch ...
    }

    /**
     * @brief Find explicit class boundary
     *
     * @param[in] k1 Class 1
     * @param[in] k2 Class 2
     *
     * @return x>value: in class 2
     */
    double findDecisionLimit(const int k1, const int k2) const {

        double min = getMean(k1);
        double max = getMean(k2);

        if (max > min)
            return
            findSingleUnivariateRootIntervalSearch(getDiscriminantFunction(k1, k2), min, max,1e-7);
        else
            return
            findSingleUnivariateRootIntervalSearch(getDiscriminantFunction(k2, k1), max, min, 1e-7);

        // @muß getestet werden
    }

    /**
     * @brief Class name for logging ...
     */
    const std::string className() const;

    /**
     * @brief Each class gets its own line, that's the easiest way
     */
    const std::string getCSVHeader() const;

};

} // ns
