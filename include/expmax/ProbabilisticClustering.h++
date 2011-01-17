/**
 * @file ProbabilisticClustering.h++
 * @author  Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 * @version 0.1
 *
 * @section LICENSE
 *
 * This file is released under LGPL v3.0.
 *
 * @section DESCRIPTION
 *
 * Gaussian Mixture Model for fitting of a known number of Gaussians;
 * and other probabilistic models of use in the processing of images.
 *
 * Common patterns are extracted afap to build a hierarchy of classes.
 *
 */

#pragma once

#include <typeinfo> // Just to print names of classes, I promise
//
//#include <boost/foreach.hpp>
//#include <boost/optional.hpp>
//
// Type assertions.
#include <boost/mpl/equal.hpp>

//#include <boost/function.hpp>
//#include <boost/lambda/bind.hpp>
//#include <boost/lambda/lambda.hpp>
//using namespace boost::lambda;

// Loop-avoiding trickery
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>

// Sine qua non of modern C++
#include <boost/shared_ptr.hpp>

//// #include <boost/array.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/symmetric.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/vector_proxy.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/ublas/lu.hpp>
//
//// Numeric limits
//#include <limits>
//
//#include <iostream> // argh, verflechtung
//
//#include "math/CircularStatistics.h++"
//#include "math/numeric_optimization.h++"
//#include "utils/comparisons.h++"
//
//#include "common_definitions.h++"

#include "GaussianMixtureModel1D.h++"

//#include "GaussianMixtureModelND.h++"
//#include "GaussianMixtureModelNDClosedForm.h++"
//#include "CircularMixtureModel1D.h++"
//
namespace CDA {

} // namespace
