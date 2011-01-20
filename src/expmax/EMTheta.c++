/**
 * @file EMTheta.c++
 * @version 0.13
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 20, 2011
 *
 */

#include "expmax/EMTheta.h++"

using namespace CDA;

EMThetas::EMThetas() {

}

EMThetas::EMThetas(const EMThetas& other) {
    std::copy(other.thetas.begin(), other.thetas.end(), std::back_inserter(thetas));
}
