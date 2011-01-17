/**
* @file GaussianMixtureModel1D.c++
*
* @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
*
*  Created on: Jan 17, 2011
*
*/

#include "expmax/ProbabilisticClustering.h++"
#include "expmax/GaussianMixtureModel1D.h++"

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "math/numeric_optimization.h++"

using namespace CDA;


inline double GaussianMixtureModel1D::squaredDistanceToMean(const unsigned k, const datapoint_t x) const {
    return fabs(x - getMean(k)); // OBS! fabs!
}

const std::string GaussianMixtureModel1D::className() const {
    return std::string("GaussianMixtureModel1D");
}

const std::string GaussianMixtureModel1D::getCSVHeader() const {
    return std::string("iter;k;p;m;sigma\n");
}

double GaussianMixtureModel1D::getMean(const unsigned k) const {
    return m_theta.getThetas()[k](1);
}

double GaussianMixtureModel1D::getSigma(const unsigned k) const {
    return m_theta.getThetas()[k](2);
}

double GaussianMixtureModel1D::evalPDF(const unsigned k, const double& x) const {
    return 1.0/(sqrt(2*M_PI) * getSigma(k)) * exp(-0.5 * squaredDistanceToMean(k, x) * pow(1 / getSigma(k), 2 ));
}

boost::function<double(const double)>
GaussianMixtureModel1D::getDiscriminantFunction(const int k1, const int k2) const {
    // @todo double-check

    assert(k1 != k2);

    return boost::function<double(const double)>(
            boost::lambda::bind(fabs, (boost::lambda::_1 - getMean(k1))) * pow(getSigma(k2),2)
            - boost::lambda::bind(fabs, (boost::lambda::_1 - getMean(k2)))  * pow(getSigma(k1),2)
            - 2 * log ( getSigma(k2) / getSigma(k1) ) * pow(getSigma(k2),2) * pow(getSigma(k1),2));
    //          < 0  ? k1 : k2 );

    // eh? ist auch falsch ...
}

double GaussianMixtureModel1D::findDecisionLimit(const int k1, const int k2) const {

    double min = getMean(k1);
    double max = getMean(k2);

    if (max > min)
        return
        findSingleUnivariateRootIntervalSearch(getDiscriminantFunction(k1, k2), min, max,1e-7);
    else
        return
        findSingleUnivariateRootIntervalSearch(getDiscriminantFunction(k2, k1), max, min, 1e-7);

}
