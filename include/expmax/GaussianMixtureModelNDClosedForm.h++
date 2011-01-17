/**
 * @file GaussianMixtureModelNDClosedForm.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 14, 2011
 *
 * @section Description
 * ... nicht wahr?
 */

#pragma once

#include "FitMultivariateMulticlassByEM.h++"

namespace CDA {

/**
 * @class GaussianMixtureModelNDClosedForm
 *
 * @brief A heteroscedastic N-dimensional Gaussian Mixture Model
 *
 */
class GaussianMixtureModelNDClosedForm : public FitMultivariateMulticlassByEM, public GaussianMixtureModelNDCommon {
    using GaussianMixtureModelNDCommon::sym_mtx_t;
    using GaussianMixtureModelNDCommon::m_cached_invsigmas;

    typedef FitMultivariateMulticlassByEM::datapoint_t datapoint_t;

public:

    /**
     * @brief Constructor
     *
     * @param[in] K_ this way round, because of default parameter in superclass
     * @param[in] D_
     */
    GaussianMixtureModelNDClosedForm(const unsigned K_, const unsigned D_)
        : FitMultivariateMulticlassByEM(K_, D_), GaussianMixtureModelNDCommon(K_, D_) {}
    // , EM<fvector_t>(D_) {}

    /**
     * @brief The quasi-MLE step.
     *
     * It is called from update_thetas
     */
    virtual void improveClusterModelParameters();

};


}
