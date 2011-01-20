/**
 * @file EM.h++
 * @version 0.12
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
 * @author flick
 *
 * @brief A generic EM model fitting class.
 *
 * @section Description
 * Contains the general outline of EM algorithm,
 * details to be filled in by implementations
 *
 * @tparam data_t The data-point holding class.
 * @tparam theta_t The model-parameter-holding class.
 */
template<class data_T, class theta_T>
class EM {

public:

    /**
     * @brief Three basic typedefs
     */
    typedef data_T data_t;
    typedef typename data_T::value_type datapoint_t;
    typedef theta_T theta_t;

protected:

    /**
     * @brief In C++, _always_ works better than inheritance. I should have remembered ...
     * Nevertheless, now you can plug in the appropriate class as a template parameter.
     */
    theta_t m_theta;

    /**
     * @brief We always need this member. Just drop it in via templating, no subclass woes
     */
    data_t m_data;

    /**
     * @brief It's here now.
     *
     * @section NOTE
     * It's virtual for historic reasons, don't replace it
     *
     * @return m_data
     */
    virtual data_t* getDataObj();

    /**
     * @brief Same in const
     *
     * @section NOTE
     * It's virtual for historic reasons, don't replace it
     *
     * @return reference to datapoints-holding object
     */
    virtual const data_t* getDataObj() const;

    /**
     * @brief It's here now.
     *
     * @section NOTE
     * It's virtual for historic reasons, don't replace it
     *
     * @return m_theta
     */
    virtual theta_t* getThetaObj();

    /**
     * @brief Same in const
     *
     * @section NOTE
     * It's virtual for historic reasons, don't replace it
     *
     * @return reference to parameter-holding object
     */
    virtual const theta_t* getThetaObj() const;

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

public:

    /**
     * @brief default constructor
     *
     * @param[in] data_t (sometimes this needs special initialization)
     * @param[in] theta_t (sometimes this needs special initialization)
     */
    EM(const data_t& data, const theta_t& theta);

    /**
     * @brief copy constructor
     * Should be implemented throughout
     */
    EM(const EM& other);

    /**
     * @brief Getter for estimate
     */
    virtual const fvector_t& getHiddenParamEstimate(const unsigned n) const = 0;

    /**
     * @brief Run Expectation Maximization
     *
     * @param[in] MAXITER up to a certain number of iterations
     * @param[in] thresh  stop when loglikelihood difference below threshold
     * @param[in] output_csv
     *
     * @todo replace the ostream with an EnhancedDataset.
     */
    void EMrun(const unsigned MAXITER,
               const double thresh,
               const boost::optional<std::ostream*> output_csv = boost::none);

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
     * @brief Initial guess for model parameters -- strictly needed
     *
     * @section NOTA
     * Input data format differs for each individual distribution,
     * be careful and read the documentation. Should be in protected
     *
     * @param[in] thetas
     */
    template<class II>
    void setTheta(std::pair<II, II> thetas) {
        // Make sure the parameters are in the right format
        BOOST_MPL_ASSERT_MSG((boost::mpl::equal<typename theta_t::value_type,
                                                typename II::value_type>::type::value),
                              UnsupportedDatavectorType, (typename II::value_type));
        m_theta . setTheta(thetas);
    }

    /**
     * @brief Implementer should overwrite this method, since typeid(this).name() returns crap
     */
    virtual const std::string className() const {
        return typeid(this).name();
    }

    /**
     * @brief Return getDataObj() . getN()
     *
     * @return Number of data points
     */
    unsigned int getN() const;

protected:

    /**
     * @brief Initialize known samples, i.e. your data points
     * Also initialize the classif vector
     */
    template<class II>
    void setData(std::pair<II, II> data_) {

        /// II must iterate over datapoint_t elements:
        // checked again in EMData
        BOOST_MPL_ASSERT_MSG((boost::mpl::equal<datapoint_t,
                              typename II::value_type>::type::value),
                              UnsupportedDatavectorType, (typename II::value_type));

#ifdef DETAIL_DEBUG_VERBOSE
        std::cerr << "EM: setData called\n";
#endif

        getDataObj() -> clear();
        getDataObj() -> setDataProper(data_);
    }

};


} // namespace
