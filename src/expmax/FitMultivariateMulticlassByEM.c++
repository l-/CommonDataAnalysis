/**
 * @file FitMultivariateMulticlassByEM.c++
 * @version 0.03
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 */

#include "expmax/FitMultivariateMulticlassByEM.h++"

using namespace CDA;

template<class theta_T>
unsigned FitMultivariateMulticlassByEM<theta_T>::getD() const {
    return getDataObj() -> getDataDimensionality();
}
