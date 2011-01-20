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
#include "VectorEMData.h++"

namespace CDA {


/**
 * @class FitMultivariateMulticlassByEM
 *
 * @brief Note that the EMData type is fixed: VectorEMData.
 *
 * @brief It's no longer a simple typedef.
 */
template<class theta_T>
class FitMultivariateMulticlassByEM : public FitMulticlassByEM<VectorEMData, theta_T> {

protected:

    /**
     * @brief We need to access this in the E(?) step
     */
    using FitMulticlassByEM<VectorEMData, theta_T>::getModifyClassif;

    /**
     * @brief Which is protected ...
     */
    using EM<VectorEMData, theta_T>::getDataObj; // which is virtual, rather that m_data

    /**
     * @same
     */
    using EM<VectorEMData, theta_T>::getThetaObj;

public:

    /**
     * @brief Define the datapoint_t
     */
    typedef VectorEMData data_t;
    typedef VectorEMData::value_type datapoint_t;
    typedef theta_T theta_t;

    /**
     * @brief Constructor
     *
     * @param[in] K_ this way round, because of default parameter in superclass
     * @param[in] data D is not explicitly passed to this class anymore since not necessary, already in data.
     * @param[in] theta
     */
    FitMultivariateMulticlassByEM(const unsigned K_, const VectorEMData& data, const theta_t& theta)
      : FitMulticlassByEM<VectorEMData, theta_t>(K_, data, theta)
    {}
    //

    unsigned getD() const;
};

} // namespace
