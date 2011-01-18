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

protected:

    EMData<datapoint_t> m_data;

public:

    /**
     * @brief Constructor
     *
     * @param[in] K_ this way round, because of default parameter in superclass
     * @param[in] D_ Data dimensionality
     */
    FitMultivariateMulticlassByEM(const unsigned K_, const unsigned D_)
    : FitMulticlassByEM<fvector_t>(K_)
    , m_data(D_)
    {}
    //

    unsigned getD() const;

    EMData<datapoint_t>* getDataObj() {
        return &m_data; }
    const EMData<datapoint_t>* getDataObj() const {
        return &m_data; }

};

} // namespace
