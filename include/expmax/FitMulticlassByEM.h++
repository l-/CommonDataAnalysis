/**
 * @file FitMulticlassByEM.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 5, 2011
 *
 */

#pragma once

#include "expmax/EM.h++"
#include <vector>

#include <boost/mpl/assert.hpp>
#include <boost/mpl/equal.hpp>

namespace CDA {

/**
 * @class FitMulticlassByEM
 *
 * An abstract model fitting class
 *
 * template param: the datapoint_t (fvector_t or double)
 *
 */
template<class T>
class FitMulticlassByEM : public EM<T> {

public:
    /**
     * @brief Propagate the datapoint_t.
     * Why can't I give it the same name? Annoying.
     */
    typedef T datapoint_t;

    using EM<datapoint_t>::className;

    /**
     * @brief We need this all over.
     */
    using EM<datapoint_t>::getN;

protected:

    /**
     * @brief Number of individual Distributions
     */
    const unsigned K;

    /**
     * @brief N classes: class appartenance probabilities
     */
    std::vector<fvector_t> classif;

    using EM<datapoint_t>::m_theta;
    using EM<datapoint_t>::getDataObj; // which is virtual, rather that m_data

    /**
     * @brief <b>E-step</b>: update hidden attributes, i.e. classification
     *
     * @section NOTE
     * It is common to all multiclass (mixture model) problem formulations,
     * thus not redefined in subclasses.
     */
    void update_hidden();

public:

    /**
     * Constructor
     *
     * To be called by subclasses
     *
     * @param[in] K_ Number of Classes (~)
     *
     */
    FitMulticlassByEM(const unsigned K_)
    : EM<datapoint_t>()
    , K(K_)
      { }

    /**
     * @brief Output CSV
     *
     * @param[in] k class
     * @param[in] iteration (optional) output iteration no. in first column
     *
     * @return one CSV line
     */
    const std::string dumpParameters(const unsigned k, const boost::optional<unsigned int> iteration) const;

    /**
     * @brief Output CSV
     *
     * @param[in] iteration (optional) output iteration no. in first column
     *
     * @return multiple CSV lines
     */
    const std::string dumpParameters(const boost::optional<unsigned int> iteration) const;

    /**
     * @brief Conveniently encapsulate access ofm_theta
     *
     * all parameters except the outer probability
     */
    double getParam(const unsigned k, const unsigned p) const;

    /**
     * @brief Header to go with the CSV output. Must be reimplemented
     */
    virtual const std::string getCSVHeader() const = 0;

    /**
     * @brief Implemented here.
     *
     * @return If you really wanna know the logLikelihood value at current iteration.
     */
    double logLikelihood() const;

    /**
     * @brief Getter for estimate
     *
     * @param[in] n Point number n
     */
    const fvector_t& getHiddenParamEstimate(const unsigned n) const;

    /**
     * @brief Better have a getter for this too
     */
    unsigned int getK() const;

    /**
     * @brief Get class with best estimated membership probability
     *
     * @param[in] n Point number n
     */
    unsigned int getBestClass(const unsigned n) const;

    /**
     * @brief Read classification results (1)
     *
     * @return Best class for each point, in the same order
     * as the points were entered (so you can just access 'em by index,
     * as everywhere)
     */
    std::vector<unsigned int> getClassifList() const;

    /**
     * @brief Get current estimated class weight
     *
     * @param[in] k class no.
     */
    double getPk(const unsigned k) const;

    /**
     * @brief Initialize known samples, i.e. your data points
     * Also initialize the classif vector
     */
    template<class II>
    void setData(std::pair<II, II> data_) {

        // II must iterate over datapoint_t elements
        BOOST_MPL_ASSERT_MSG((boost::mpl::equal<datapoint_t, typename II::value_type>::type::value), UnsupportedDatavectorType, (typename II::value_type));

#ifdef VERBOSE
        std::cerr << "FitMulticlassByEM: setData called\n";
#endif

        getDataObj() . setDataProper(data_);
        initClassif();
    }

    /**
     * @brief Initialize class membership beliefs
     *
     * @section NOTE
     * setData of derived classes must call this
     * => done, defined setData here. can one be "using" a template fn?
     */
    void initClassif();

    /**
     * @brief Evaluate the cluster model PDF with class parameters \f$k\f$ at \f$\vec{x}\f$
     *
     * @param[in] k class
     * @param[in] x feature vector
     *
     * @return probability \f$p(x,k\vert\theta)\f$
     */
    virtual double evalPDF(const unsigned k, const datapoint_t& x) const = 0;

};


} // namespace
