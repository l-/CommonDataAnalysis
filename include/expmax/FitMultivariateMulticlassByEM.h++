/**
 * @file FitMultivariateMulticlassByEM.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 5, 2011
 *
 */

#pragma once

#include "FitMulticlassByEM.h++"

namespace CDA {

/**
 * @class FitMultivariateMulticlassByEM
 *
 * @brief It's no longer a simple typedef.
 */
class FitMultivariateMulticlassByEM : public FitMulticlassByEM<fvector_t> {

protected:
    using FitMulticlassByEM<fvector_t>::classif;
    using FitMulticlassByEM<fvector_t>::getN;

public:

    /**
     * @brief Define the datapoint_t
     */
    typedef fvector_t datapoint_t;

    /**
     * @brief get number of classes
     */
    using FitMulticlassByEM<fvector_t>::getK;

public:

    /**
     * @brief Constructor
     *
     * @param[in] K_ this way round, because of default parameter in superclass
     * @param[in] D_
     */
    FitMultivariateMulticlassByEM(const unsigned K_, const unsigned D_)
    : FitMulticlassByEM<fvector_t>(K_, D_) {}
    // , EM<fvector_t>(D_) {}

    /**
     * @brief Convenience getter for D, same format as getN() and getP()
     *
     * @return Data space dimensionality
     */
    unsigned getD() const;
};

} // namespace
