/**
 * @file GaussianMixtureModelNDCommon.c++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 */

#include "expmax/GaussianMixtureModelNDCommon.h++"

#include "common_definitions.h++"

#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp> // for vector_range
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include "math/matrix_functions.h++"

#include <cmath> // hallo?

using namespace CDA;

GaussianMixtureModelNDCommon::
GaussianMixtureModelNDCommon(const unsigned K_, const unsigned D_, const data_t& data, const theta_t& theta)
    // : EMData<fvector_t>(D_) // <-- don't forget!!!
    : EMGenericMixtureModelCore<GaussianMixtureModelNDParams>(K_, 1+D_+(D_*(D_+1))/2, D_, data, theta)
{
#ifdef DETAIL_VERBOSE_2
        std::cerr << "GaussianMixtureModelNDCommon: constructor was called. P=" << getP() << " " << getK() << " " << getD() << std::endl;
#endif
}

double GaussianMixtureModelNDCommon::evalPDF(const unsigned k, const fvector_t& x) const {

    assert(m_theta . getThetas() . size() > 0);

    namespace ublas = boost::numeric::ublas;

    ublas::vector<double> zwe(getD()); // vector Zwischenergebnis
    ublas::vector<double> xmu = x - getThetaObj() -> getMean(k);
    zwe = ublas::prod(getThetaObj() -> getCachedInvSigma(k), xmu);

    double res = pow(2*M_PI, -0.5*(double)getD()) *
                    pow(getThetaObj() -> getCachedSigmaDet(k), -0.5) *
                    exp(- 0.5 *  ublas::inner_prod(xmu, zwe));

#ifdef DETAIL_VERBOSE_2
    if (isnan(res)) {
        std::cout << "Asked to evalPDF at sigma " << getThetaObj() -> getSigmaMatrix(k) << ", mean " << getThetaObj() -> getMean(k) << std::endl;
        std::cout << "(1) " << pow(2*M_PI, -0.5*(double)getD()) << " " << pow(getThetaObj() -> getCachedSigmaDet(k), -0.5) << std::endl;
        std::cout << "(1) " << xmu << std::endl;
        std::cout << "(1) " << getThetaObj() -> getCachedSigmaDet(k) << " " << getD() << std::endl;
        std::cout << "(1) " << zwe << std::endl;
        std::cout << "(2) " << exp(- 0.5 *  ublas::inner_prod(xmu, zwe)) << std::endl;
        std::cout << pow(2*M_PI, -0.5*(double)getD()) * pow(getThetaObj() -> getCachedSigmaDet(k), -0.5) * exp(- 0.5 *  ublas::inner_prod(xmu, zwe)) << std::endl;
    }
#endif

    return res;
}

size_t GaussianMixtureModelNDCommon::getD() const {
    return getDataObj() -> getDataDimensionality();
}
