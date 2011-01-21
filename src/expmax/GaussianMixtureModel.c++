/**
 * @file GaussianMixtureModel.c++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 * n-dimensional Gaussian Mixture MOdel. common "data" core minus the procedures for optimization
 *
 */

#include "expmax/GaussianMixtureModelNDCommon.h++"
#include "expmax/GaussianMixtureModel.h++"

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>

using namespace CDA;

const std::string GaussianMixtureModel::paramName(const unsigned p) const {
    return GaussianMixtureModelNDCommon::paramName(p);
}

const std::string GaussianMixtureModel::className() const {
    return std::string("GaussianMixtureModel(NDClosedForm)");
}

void GaussianMixtureModel::improveClusterModelParameters() {

#ifdef VERBOSE
    std::cout << className() << ": M step.\n";
#endif

    unsigned K = getK();

    double sumclassif[K];
    for (unsigned k=0; k<K; ++k) {
        sumclassif[k] = 0.0;
    }

    // Optimize all you like, someday

    std::vector<fvector_t> means;
    std::vector<sym_mtx_t> sigmas;
    for (unsigned k=0; k<K; ++k) {
        means.push_back(boost::numeric::ublas::zero_vector<double>(getD()));
        sigmas.push_back(
                sym_mtx_t(
                        boost::numeric::ublas::zero_matrix<double>(getD(), getD())));
    }

    // Preparation
    for (unsigned n=0; n<getN(); ++n) {
        for (unsigned k=0; k<getK(); ++k) {
            sumclassif[k] += getClassif(n)(k);
            means[k] += getClassif(n)(k) * getDataObj() -> getData(n);
        }
    }

    // Calculate weighted sample mean, in effect
    for (unsigned k=0; k<K; ++k) {

#ifdef DETAIL_VERBOSE_2
    std::cout << "Sumclassif " << k << " = " << sumclassif[k] << " .\n";
#endif


        means[k] /= sumclassif[k];

        getThetaObj() -> getModifyMean(k) = means[k];
//
//        for (unsigned d=0; d<getD(); ++d) {
//
//            m_theta.getModifyThetas()[k][1+d] = means[k](d);
//
#ifdef DETAIL_VERBOSE_2
        std::cout << "Mean " << k << " now " << getThetaObj() -> getMean(k) << std::endl;
#endif
//
//        }
    }

    // Calculate weighted sample covariance matrix, in effect
    for (unsigned n=0; n<getN(); ++n) {
        for (unsigned d=0; d<getD(); ++d) {
            for (unsigned e=d; e<getD(); ++e) {
                for (unsigned k=0; k<getK(); ++k) {
                    sigmas[k](d,e) += classif[n](k) *
                           ((getDataObj()->getData(n)(d) - means[k](d)) * (getDataObj()->getData(n)(e) - means[k](e)));
                }
            }
        }
    }

    // Somehow, I doubt this is correct ...
    // getThetaObj() ->



    for (unsigned k=0; k<getK(); ++k) {
        setSigma(k, sigmas[k]/sumclassif[k]);
#ifdef DETAIL_VERBOSE_2
        std::cout << "Sigmas" << " " << k << " now " << sigmas[k]/sumclassif[k] << std::endl;
#endif
    }


    getThetaObj() -> updateCached();

    // @todo: improve efficiency and calculations.

}

