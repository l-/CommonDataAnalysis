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
#include "math/matrix_functions.h++"

using namespace CDA;

GaussianMixtureModelNDCommon::
GaussianMixtureModelNDCommon(const unsigned K_, const unsigned D_)
    : EMData<fvector_t>(D_) // <-- don't forget!!!
    , EMGenericMixtureModelCore(K_, D_, D_+D_*(D_+1)/2)

{

}

const std::string GaussianMixtureModelNDCommon::paramName(const unsigned p) const {
    std::stringstream out;
    if (p < getDataDimensionality()) {
        out << "m_" << p;
    } else {
        int a = p - getDataDimensionality();
        out << "s_" << i(a) << "_" << j(a);
    }
    return out.str();
}


const boost::numeric::ublas::vector_range<const fvector_t>
GaussianMixtureModelNDCommon::getMean(const unsigned k) const {
    using namespace boost::numeric::ublas;
    return project(m_theta.getThetas()[k], range(1, 1+getD()));
}

EMData<GaussianMixtureModelNDCommon::datapoint_t>& GaussianMixtureModelNDCommon::getDataObj() {
    // Let's hope it worx
    return *this;
}

const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper>
GaussianMixtureModelNDCommon::getInvSigma(const unsigned k) const {

    // @todo "calculateInvSigma"

    namespace ublas = boost::numeric::ublas;
    // @todo avoid recalculations. but first get it right

    ublas::matrix<double> inverse(getD(), getD());
    typedef ublas::permutation_matrix<std::size_t> permatrix;

    // create a working copy of the input, as a matrix (can't be symmetric_matrix, na?)
    ublas::matrix<double> A = getSigmaMatrix(k);

    // create a permutation matrix for the LU-factorization
    permatrix pm(A.size1());

    // perform LU-factorization with pivoting
    lu_factorize(A,pm);

    // create identity matrix of "inverse"
    inverse.assign(ublas::identity_matrix<double>(A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);
    return inverse;
}

double GaussianMixtureModelNDCommon::getSigmaDet(const unsigned k) const {
    return det(getSigmaMatrix(k));
}

const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper>
GaussianMixtureModelNDCommon::getSigmaMatrix(const unsigned k) const {
    boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> Sigma(getD(), getD());
    // @todo c'mon, just copy the memory
    for (unsigned p = getD(); p < P; ++p) {
        unsigned a = p - getD();
        unsigned i_ = i(a);
        unsigned j_ = j(a);
        Sigma(i_,j_) = getParam(k, p);
    }
    return Sigma;
}

double GaussianMixtureModelNDCommon::evalPDF(const unsigned k, const fvector_t& x) const {

    assert(m_theta . getThetas() . size() > 0);

    namespace ublas = boost::numeric::ublas;

    ublas::vector<double> zwe(getD()); // vector Zwischenergebnis
    ublas::vector<double> xmu = x - getMean(k);
    zwe = ublas::prod(getInvSigma(k), xmu);

//    std::cout << "Asked to evalPDF at sigma " << getSigmaMatrix(k) << ", mean " << getMean(k) << std::endl;
//    std::cout << "(1) " << pow(2*M_PI, -0.5*(double)getD()) * 1/sqrt(getSigmaDet(k)) << std::endl;
//    std::cout << "(1) " << xmu << std::endl;
//    std::cout << "(1) " << getSigmaDet(k) << std::endl;
//    std::cout << "(1) " << zwe << std::endl;
//    std::cout << "(2) " << exp(- 0.5 *  ublas::inner_prod(xmu, zwe)) << std::endl;
//    std::cout << pow(2*M_PI, -0.5*(double)getD()) * 1/sqrt(getSigmaDet(k)) * exp(- 0.5 *  ublas::inner_prod(xmu, zwe)) << std::endl;

    return pow(2*M_PI, -0.5*(double)getD()) * 1/sqrt(getSigmaDet(k)) * exp(- 0.5 *  ublas::inner_prod(xmu, zwe));

    // Seems OK
}
