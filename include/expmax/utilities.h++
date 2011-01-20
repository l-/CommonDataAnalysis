/**
 * @file utilities.h++
 * @version 0.14
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 * Some commonly used utility functions based on the ExpMax classes
 */

#pragma once

#include "expmax/ProbabilisticClustering.h++"

#include <boost/algorithm/minmax_element.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace CDA {

/* UTILITIES: */

/**
 * @brief
 * Use univariate mixture model to determine a threshold
 * The classical application, from BV1 lecture actually!
 */
template<class Iter>
boost::shared_ptr<GaussianMixtureModel1D> thresholdFinder(const std::pair<Iter, Iter> input,
        double bgprob = 0.5,
        const boost::optional<std::ostream*> output_csv = boost::none
) {
    std::vector<fvector_t> init_theta_bimodal;

    fvector_t hintergrund(3);
    fvector_t vordergrund(3);

    std::pair<Iter, Iter> minmax =
            boost::minmax_element(input.first, input.second);

    // Init prob.
    hintergrund(0) = bgprob; // if you know differently, please initialize it accordingly
    vordergrund(0) = 1 - hintergrund(0);

    // Init mean
    hintergrund(1) = * ( minmax.first );
    vordergrund(1) = * ( minmax.second );

    // Init stddev
    hintergrund(2) = 0.5 * (*minmax.second - *minmax.first);
    vordergrund(2) = hintergrund(2);

    init_theta_bimodal.push_back(hintergrund);
    init_theta_bimodal.push_back(vordergrund);

    GaussianMixtureModel1D::data_t data;
    data.setDataProper(input);
    GaussianMixtureModel1D::theta_t theta;
    theta.setTheta(std::make_pair(init_theta_bimodal.begin(), init_theta_bimodal.end()));

    boost::shared_ptr<GaussianMixtureModel1D> bimodal (new GaussianMixtureModel1D(2, data, theta ) );

    bimodal -> EMrun(50, 0.0000001, output_csv); // @todo time

    return bimodal;
}

} // namespace
