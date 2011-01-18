/**
 * @file FitUnivariateMulticlassByEM.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 */

#pragma once

#include "expmax/FitMulticlassByEM.h++"

namespace CDA {

/**
 * @brief It's been demoted to a typedef.
 */
class FitUnivariateMulticlassByEM : public FitMulticlassByEM<double> {

    EMData<double> m_data;

public:

    /**
     * @brief Constructor
     *
     * @param[in] K_ this way round, because of default parameter in superclass
     */
    FitUnivariateMulticlassByEM(const unsigned K_)
    : FitMulticlassByEM<double>(K_)
    , m_data(1)
    {}
    //
};

} // namespace
