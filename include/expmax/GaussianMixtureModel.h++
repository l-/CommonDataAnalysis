/**
 * @file GaussianMixtureModel.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 14, 2011
 *
 * @section Description
 * ... nicht wahr?
 */

#pragma once

#include "expmax/GaussianMixtureModelNDCommon.h++"
#include "expmax/FitMultivariateMulticlassByEM.h++"

namespace CDA {

/**
 * @class GaussianMixtureModelNDClosedForm
 *
 * @brief A heteroscedastic N-dimensional Gaussian Mixture Model
 *
 */
class GaussianMixtureModel: public GaussianMixtureModelNDCommon {

//    using GaussianMixtureModelNDCommon::getClassif;
//
//    using GaussianMixtureModelNDCommon::getP; // zufrieden?

public:

    /**
     * @brief Types
     */
    typedef GaussianMixtureModelNDCommon::datapoint_t datapoint_t;
    /**
     * @brief Types
     */
    typedef GaussianMixtureModelNDCommon::data_t data_t;
    /**
     * @brief Types
     */
    typedef GaussianMixtureModelNDCommon::theta_t theta_t;
    /**
     * @brief Types
     */
    typedef GaussianMixtureModelNDCommon::theta_t::sym_mtx_t sym_mtx_t;

    /**
     * @brief Constructor
     *
     * @param[in] K_ this way round, because of default parameter in superclass
     * @param[in] D_
     * @param[in] data
     * @param[in] theta
     */
    GaussianMixtureModel(const unsigned K_, const unsigned D_, const data_t& data, const theta_t& theta)
        : GaussianMixtureModelNDCommon(K_, D_, data, theta) {}

    /**
     * @brief Constructor, initializing data etc. appropriately
     *
     * @param[in] K_ this way round, because of default parameter in superclass
     * @param[in] D_
     *
     */
    GaussianMixtureModel(const unsigned K_, const unsigned D_)
            : GaussianMixtureModelNDCommon(K_, D_, data_t(D_), theta_t(K_, D_)) {}

//    template<class II>
//    void initializeClusters(std::pair<II, II> thetas_) {
//        getThetaObj() -> setTheta(thetas_);
//        getDataObj() -> updateCached();
//    }

    // this is turning ugly once more and not what I want when I order inheritance
    using GaussianMixtureModelNDCommon::getSigma;
    using GaussianMixtureModelNDCommon::getCachedSigmaDet;
    using GaussianMixtureModelNDCommon::getMean;
    using GaussianMixtureModelNDCommon::getClassProb;

protected:

    /**
     * @brief The quasi-MLE step, part of <b>E-step</b>.
     *
     * @section CALLED-FROM
     * It is called from update_thetas, which is already defined
     * in EMGenericMixtureModelCore.c++
     */
    virtual void improveClusterModelParameters();

    /**
     * @brief Construct names of parameters, e.g. for output
     * calls same function from ...Common
     */
    virtual const std::string paramName(const unsigned p) const;

    /**
     * @brief Class name for logging
     * @return
     */
    virtual const std::string className() const;

};


}
