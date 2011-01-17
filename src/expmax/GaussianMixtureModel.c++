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
    std::vector<boost::numeric::ublas::symmetric_matrix<double> > sigmas;
    for (unsigned k=0; k<K; ++k) {
        means.push_back(boost::numeric::ublas::zero_vector<double>(getDataDimensionality()));
        sigmas.push_back(boost::numeric::ublas::zero_matrix<double>(getDataDimensionality(), getDataDimensionality()));
    }

    // TODO

    for (unsigned n=0; n<getN(); ++n) {
        for (unsigned k=0; k<getK(); ++k) {
            sumclassif[k] += classif[n](k);
            means[k] += classif[n](k) * getData(n);
        }
    }

    for (unsigned k=0; k<K; ++k) {
        means[k] /= sumclassif[k];

        for (unsigned d=0; d<getDataDimensionality(); ++d) {

            m_theta.getModifyThetas()[k][1+d] = means[k](d);

#ifdef VERBOSE
        std::cout << "Mean_" << d << " " << k << " now " << m_theta.getThetas()[k](0) << std::endl;
#endif


        }
    }


//    for (unsigned n=0; n<m_data.getN(); ++n) {
//        for (unsigned k=0; k<K; ++k) {
//            // using NEW mean, thus in extra step.
//
//            // @todo: effizienter und numerisch besser berechnen ... aber erstmal das ganze verfahren geradebiegen.
//            m_theta.getModifyThetas()[k](2) +=
//                classif[n](k) * squaredDistanceToMean(k, m_data.getData(n));
//        }
//    }
//
//    for (unsigned k=0; k<K; ++k) {
//        m_theta.getModifyThetas()[k](2) = sqrt(m_theta.getModifyThetas()[k](2) / sumclassif[k]);
//        m_theta.getModifyThetas()[k](0) = 1/((double)(m_data.getN())) * sumclassif[k];
//    }

    // @todo atomic
}

