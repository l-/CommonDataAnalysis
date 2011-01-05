/**
 * @file numeric_optimization.h++
 *
 * @version 0.01
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Dec 28, 2010
 *
 */

#pragma once

#include <boost/function.hpp>

namespace CDA {

/**
 * @function findSingleUnivariateRootIntervalSearch
 *
 * @brief Finds one root in the given interval
 * @param[in] f the function object
 * @param[xa] initial lower bound for x
 * @param[xb] initial upper bound for x
 * @param[epsilon] get this near, x-wise
 *
 * @return some value with f(y) = x, hopefully. (also if there's a plateau).
 * or a value with a discontinuity and different signs on both sides
 *
 */
double findSingleUnivariateRootIntervalSearch(const boost::function<double(double)>& f,
                                              const double xa,
                                              const double xb,
                                              const double epsilon);

};
