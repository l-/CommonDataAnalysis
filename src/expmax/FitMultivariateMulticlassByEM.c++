/**
 * @file FitMultivariateMulticlassByEM.c++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 */

#include "expmax/FitMultivariateMulticlassByEM.h++"

using namespace CDA;

unsigned FitMultivariateMulticlassByEM::getD() const {
    return getDataObj() . getDataDimensionality();
}
