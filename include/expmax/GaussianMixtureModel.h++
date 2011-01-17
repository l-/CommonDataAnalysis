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

#include "GaussianMixtureModelNDCommon.h++"
#include "FitMultivariateMulticlassByEM.h++"

namespace CDA {

/**
 * @class GaussianMixtureModelNDClosedForm
 *
 * @brief A heteroscedastic N-dimensional Gaussian Mixture Model
 *
 */
class GaussianMixtureModel: public GaussianMixtureModelNDCommon {
    // Fit ...
    using GaussianMixtureModelNDCommon::sym_mtx_t;
    using GaussianMixtureModelNDCommon::m_cached_invsigmas;

    typedef GaussianMixtureModelNDCommon::datapoint_t datapoint_t;

    using GaussianMixtureModelNDCommon::classif;

public:

    /**
     * @brief Constructor
     *
     * @param[in] K_ this way round, because of default parameter in superclass
     * @param[in] D_
     */
    GaussianMixtureModel(const unsigned K_, const unsigned D_)
        : GaussianMixtureModelNDCommon(K_, D_) {}
    // FitMultivariateMulticlassByEM(K_, D_),
    // , EM<fvector_t>(D_) {}

protected:

    /**
     * @brief The quasi-MLE step, part of <b>E-step</b>.
     *
     * It is called from update_thetas, which is already defined
     * in EMGenericMixtureModelCore.c++
     */
    virtual void improveClusterModelParameters();

    /**
     * @brief Construct names of parameters, e.g. for output
     * calls same function from ...Common
     */
    virtual const std::string paramName(const unsigned p) const;

};


}
