/**
 * @file EM.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 5, 2011
 *
 */

#pragma once

#include "EMData.h++"
#include "EMTheta.h++"
#include "common_definitions.h++"

namespace CDA {

/**
 * @class EM
 *
 * @brief A generic EM model fitting class.
 * Contains the general outline of EM algorithm,
 * details to be filled in by implementations (not that easy in C++ ...)
 */
template<class datapoint_t>
class EM {

protected:

    /**
     * @brief In C++, _always_ works better than inheritance. I should have remembered ...
     */
    EMThetas m_theta;

    /**
     * @brief In C++, _always_ works better than inheritance
     */
    EMData<datapoint_t> m_data;

    /**
     * @brief One iteration.
     */
    double EMstep() {
        update_hidden();
        update_thetas();
        return logLikelihood();
    }

    /**
     * @brief <b>E-step</b>: update hidden attributes, i.e. class membership probabilities
     */
    virtual void update_hidden() = 0;

    /**
     * @brief <b>M-step</b>: update parameters of model PDF
     */
    virtual void update_thetas() = 0;

    /**
     * @brief Important. This is being optimized
     */
    virtual double logLikelihood() const = 0;

    /**
     * @brief Implementer should overwrite this method, since typeid(this).name() returns crap
     */
    virtual const std::string className() const {
        return typeid(this).name();
    }

public:

    /**
     * @brief Constructor, to be called explicitly in multivariate case
     */
    EM(const unsigned D_ = 1)
    : m_data(D_) { }

    /**
     * @brief Getter for estimate
     */
    virtual const fvector_t& getHiddenParamEstimate(const unsigned n) const = 0;

    /**
     * @brief Initial guess for model parameters -- strictly needed
     *
     * Input data format differs for each individual distribution,
     * be careful and read the documentation.
     *
     * @param[in] thetas
     */
    template<class II>
    void setTheta(std::pair<II, II> thetas) {
        m_theta . setTheta(thetas);
    }


    /**
     * @brief Run Expectation Maximization
     *
     * @param[in] MAXITER up to a certain number of iterations
     * @param[in] thresh  stop when loglikelihood difference below threshold
     * @param[in] output_csv
     *
     */
    void EMrun(const unsigned MAXITER, const double thresh, const boost::optional<std::ostream*> output_csv = boost::none);

    /**
     * @brief Output CSV
     *
     * @param[in] iteration (optional) output iteration no. in first column
     *
     */
    virtual const std::string dumpParameters(const boost::optional<unsigned int> iteration) const = 0;

    /**
     * @brief Möchte es nicht zu kompliziert machen (14/12/2010),
     * daher bekommen die ganzen Klassen erst einmal CSV-Ausgabe hardwired.
     */
    virtual const std::string getCSVHeader() const = 0;

    /**
     * @brief Initialize known samples, i.e. your data points
     * Also initialize the classif vector
     */
    template<class II>
    void setData(std::pair<II, II> data_) {

        // II must iterate over datapoint_t elements
        BOOST_MPL_ASSERT_MSG((boost::mpl::equal<datapoint_t, typename II::value_type>::type::value), UnsupportedDatavectorType, (typename II::value_type));

        m_data . setDataProper(data_);
    }
};



} //ns
