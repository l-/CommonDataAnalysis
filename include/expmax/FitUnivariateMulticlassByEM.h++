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
class FitUnivariateMulticlassByEM : public FitMulticlassByEM<EMData<double>, EMThetas> {

public:

    typedef double datapoint_t;
    typedef EMData<double> data_t;
    typedef EMThetas theta_t;

    /**
     * @brief Constructor
     *
     * @param[in] K_ this way round, because of default parameter in superclass
     * @param[in] data
     * @param[in] theta
     */
    FitUnivariateMulticlassByEM(const unsigned K_, const data_t& data, const theta_t& theta)
      : FitMulticlassByEM<EMData<double>, EMThetas>(K_, data, theta)
    // , m_data(1)
    {}
    //
};

} // namespace
