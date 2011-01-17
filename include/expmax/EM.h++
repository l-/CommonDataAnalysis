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

#include <boost/optional.hpp>

// Type assertions to check input types
#include <boost/mpl/equal.hpp>
#include <boost/mpl/assert.hpp>

namespace CDA {

/**
 * @class EM
 *
 * @brief A generic EM model fitting class.
 *
 * @section Description
 * Contains the general outline of EM algorithm,
 * details to be filled in by implementations
 */
template<class datapoint_t>
class EM {

protected:

    /**
     * @brief In C++, _always_ works better than inheritance. I should have remembered ...
     */
    EMThetas m_theta;

    /**
     * @brief In C++, has-a works better than is-a ...
     *
     * @section NOTA
     * Careful, only access via accessor
     */
    EMData<datapoint_t> m_data;

    /**
     * @brief Might be redefined by a subclass in case it wants its own version of EMData ...
     * so be careful, ONLY access it via accessor
     *
     * @return m_data
     */
    virtual EMData<datapoint_t>& getDataObj();

    /**
     * @brief Same in const
     * @return reference to datapoints-holding object
     */
    virtual const EMData<datapoint_t>& getDataObj() const;

    /**
     * @brief Return getDataObj() . getN()
     *
     * @return Number of data points
     */
    unsigned int getN() const;

    /**
     * @brief One iteration.
     *
     * @section DESCRIPTION
     * This is the basic EM (not extended EM etc.) algorithm,
     * by virtue of calls to virtual functions you can adapt
     * the procedure to your problem.
     *
     * @section SIDE-EFFECTS
     * This will update the model parameters stored in the EMThetas class.
     *
     * @return The new log-likelihood value after performing the iteration.
     */
    double EMstep() {
        update_hidden();
        update_thetas();
        return logLikelihood();
    }

    /**
     * @brief <b>E-step</b>: update hidden attributes
     */
    virtual void update_hidden() = 0;

    /**
     * @brief <b>M-step</b>: update accessible parameters of model PDF
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
     * @brief Constructor, to be called explicitly in multivariate case only
     */
    EM(const unsigned D_ = 1);

    /**
     * @brief Getter for estimate
     */
    virtual const fvector_t& getHiddenParamEstimate(const unsigned n) const = 0;

    /**
     * @brief Initial guess for model parameters -- strictly needed
     *
     * @section NOTA
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

        /// II must iterate over datapoint_t elements:
        BOOST_MPL_ASSERT_MSG((boost::mpl::equal<datapoint_t, typename II::value_type>::type::value), UnsupportedDatavectorType, (typename II::value_type));

        m_data . setDataProper(data_);
    }
};


} // namespace
