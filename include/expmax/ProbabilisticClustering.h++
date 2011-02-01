/**
 * @file ProbabilisticClustering.h++
 * @author  Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 * @version 0.1
 *
 * @section LICENSE
 *
 * see note ...
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

// Just to print names of classes, I promise
#include <typeinfo>

// Type assertions.
#include <boost/mpl/equal.hpp>

// Loop-avoiding trickery
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>

// Sine qua non of modern C++
#include <boost/shared_ptr.hpp>

#include "GaussianMixtureModel1D.h++"
#include "GaussianMixtureModel.h++"

namespace CDA {

} // namespace
