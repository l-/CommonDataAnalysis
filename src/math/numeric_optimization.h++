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

/**
 * @function findSingleUnivariateRootIntervalSearch
 *
 * Finds one root
 *
 */
double findSingleUnivariateRootIntervalSearch(const boost::function<double(double)>& f,
                                              const double xa,
                                              const double xb,
                                              const double epsilon);
