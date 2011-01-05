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

#include <boost/foreach.hpp>
#include <boost/optional.hpp>

// Type assertions.
#include <boost/mpl/equal.hpp>

#include <boost/function.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
using namespace boost::lambda;

// Loop-avoiding trickery
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/algorithm/minmax_element.hpp>

// Sine qua non of modern C++
#include <boost/shared_ptr.hpp>

// #include <boost/array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

// Numeric limits
#include <limits>

#include <iostream> // argh, verflechtung

#include "math/CircularStatistics.h++"
#include "math/numeric_optimization.h++"
#include "utils/comparisons.h++"

#include "common_definitions.h++"
#include "GaussianMixtureModel1D.h++"
#include "GaussianMixtureModelND.h++"
#include "CircularMixtureModel1D.h++"

namespace CDA {

/* UTILITIES: */

/**
 * @brief
 * Use univariate mixture model to determine a threshold
 * The classical application, from BV1 lecture actually!
 */
template<class Iter>
boost::shared_ptr<GaussianMixtureModel1D> thresholdFinder(const std::pair<Iter, Iter> input,
        const boost::optional<std::ostream*> output_csv = boost::none
) {

    std::vector<fvector_t> init_theta_bimodal;

    fvector_t hintergrund(3);
    fvector_t vordergrund(3);

    std::pair<Iter, Iter> minmax =
            boost::minmax_element(input.first, input.second);

    // Init prob.
    hintergrund(0) = 0.5; // if you know differently, please initialize it accordingly
    vordergrund(0) = 1 - hintergrund(0);

    // Init mean
    hintergrund(1) = * ( minmax.first );
    vordergrund(1) = * ( minmax.second );

    // Init stddev
    hintergrund(2) = 0.5 * (*minmax.second - *minmax.first);
    vordergrund(2) = hintergrund(2);

    init_theta_bimodal.push_back(hintergrund);
    init_theta_bimodal.push_back(vordergrund);

    boost::shared_ptr<GaussianMixtureModel1D> bimodal (new GaussianMixtureModel1D(2) );

    bimodal -> setData(input);
    bimodal -> setTheta(std::make_pair(init_theta_bimodal.begin(), init_theta_bimodal.end()));
    bimodal -> EMrun(50, 0.0000001, output_csv); // @todo time

    return bimodal;
}


/**
 * @brief Calculate the determinant of a small matrix.
 */
template<class matrix_t>
double det(matrix_t& m) {

    if (m.size1() == 1) {
        return m(0,0);
    }

    namespace ublas = boost::numeric::ublas;

    typedef ublas::permutation_matrix<std::size_t> permatrix;

    ublas::matrix<double> A = m;
    permatrix pm(A.size1());
    int lures = ublas::lu_factorize(A, pm);

    if( lures != 0 ) {
        std::cerr << "LU factorization of " << m << " failed." << std::endl;
        std::cerr << "you are stuck with: " << A << pm << std::endl;
    }

    double res = 1.0;
    for (unsigned i=0; i<A.size1(); ++i) {
        res *= A(i,i);
    }
    return res;
}


///**
// * Note: Circular and other manifold statistics: Maybe faking them via Gausses
// * (and specially defined distances) is en
// */

} // namespace
