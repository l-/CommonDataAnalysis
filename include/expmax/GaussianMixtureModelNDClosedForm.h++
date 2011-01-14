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

};


}
