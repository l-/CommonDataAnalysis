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
            sumclassif[k] += classif[n](k);
            means[k] += classif[n](k) * getData(n);
        }
    }

    // Calculate weighted sample mean, in effect
    for (unsigned k=0; k<K; ++k) {
        means[k] /= sumclassif[k];

        for (unsigned d=0; d<getD(); ++d) {

            m_theta.getModifyThetas()[k][1+d] = means[k](d);

#ifdef VERBOSE
        std::cout << "Mean_" << d << " " << k << " now " << m_theta.getThetas()[k](1+d) << std::endl;
#endif

        }
    }

    // Calculate weighted sample covariance matrix, in effect
    for (unsigned n=0; n<getN(); ++n) {
        for (unsigned d=0; d<getD(); ++d) {
            for (unsigned e=d; e<getD(); ++e) {
                for (unsigned k=0; k<getK(); ++k) {
                    sigmas[k](d,e) += classif[n](k) *
                           ((getData(n)(d) - means[k](d)) * (getData(n)(e) - means[k](e)));
                }
            }
        }
    }

    // Somehow, I doubt this is correct ...
    for (unsigned k=0; k<getK(); ++k) {
        for (unsigned d=0; d<getDataDimensionality(); ++d) {
            for (unsigned e=d; e<getDataDimensionality(); ++e) {

#ifdef EXTRA_VERBOSE
                std::cerr << d << " " << e << " " << 1 + getD() + a(e,d) << std::endl;
#endif

                m_theta.getModifyThetas()[k][1 + getD() + a(e,d)] = sigmas[k](d,e) / sumclassif[k];

            }
        }

    #ifdef VERBOSE
            std::cout << "Sigmas" << " " << k << " now " << sigmas[k]/sumclassif[k] << std::endl;
    #endif
    }



    // @todo: effizienter und numerisch besser berechnen ... aber erstmal das ganze verfahren geradebiegen.

    // @todo atomic
}

