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
GaussianMixtureModelNDCommon(const unsigned K_, const unsigned D_)
    // : EMData<fvector_t>(D_) // <-- don't forget!!!
    : EMGenericMixtureModelCore(K_, 1+D_+(D_*(D_+1))/2, D_)
{

#ifdef DETAIL_VERBOSE_2
        std::cerr << "GaussianMixtureModelNDCommon: constructor was called. P=" << getP() << " " << getK() << " " << getD() << std::endl;
#endif
}

const std::string GaussianMixtureModelNDCommon::paramName(const unsigned p) const {
    std::stringstream out;
    if (p < getD()) {
        out << "m_" << p;
    } else {
        int a = p - getD(); // getDataDimensionality should work just as well!!!
        out << "s_" << i(a) << "_" << j(a);
    }
    return out.str();
}


const boost::numeric::ublas::vector_range<const fvector_t>
GaussianMixtureModelNDCommon::getMean(const unsigned k) const {
    using namespace boost::numeric::ublas;
    return project(m_theta.getThetas()[k], range(1, 1+getD()));
}

const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper>
GaussianMixtureModelNDCommon::getCachedInvSigma(const unsigned k) const {
    return m_cached_invsigmas[k];
}

double
GaussianMixtureModelNDCommon::getCachedSigmaDet(const unsigned k) const {
    return m_cached_sigmadet[k];
}

void GaussianMixtureModelNDCommon::updateCached() {
    m_cached_invsigmas.clear();
    m_cached_sigmadet.clear();
    for (unsigned k=0; k<getK(); ++k) {
        m_cached_invsigmas.push_back(getInvSigma(k));
        m_cached_sigmadet.push_back(getSigmaDet(k));
        // Now only once per turn
#ifdef VERBOSE
        std::cerr << "SIGMA " << getSigmaMatrix(k) << std::endl;
        std::cerr << "INV SIGMA " << getInvSigma(k) << std::endl;
        std::cerr << "DET(SIGMA)" << getSigmaDet(k) << std::endl;
#endif
    }
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

    try {
        try {
            lu_substitute(A, pm, inverse);
        } catch ( boost::numeric::ublas::internal_logic& e ) {
            std::cerr << e.what() << std::endl;
        }
    } catch ( std::logic_error& e ) {
        std::cerr << e.what() << "\n";
    }

    // Zufrieden?

    sym_mtx_t result;
    try {
        result = boost::numeric::ublas::triangular_adaptor<boost::numeric::ublas::matrix<double>, boost::numeric::ublas::upper>(inverse);
        // Is there a problem? Just give ZERO
    } catch ( boost::numeric::ublas::external_logic& e ) {
        std::cerr << e.what() << std::endl;
        // result = boost::numeric::ublas::identity_matrix<double>(getD());
        result = boost::numeric::ublas::zero_matrix<double>(getD(), getD());
    }

    return result;
}

double GaussianMixtureModelNDCommon::getSigmaDet(const unsigned k) const {
    return det(getSigmaMatrix(k));
}

const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper>
GaussianMixtureModelNDCommon::getSigmaMatrix(const unsigned k) const {
    boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> Sigma(getD(), getD());
    // @todo c'mon, just copy the memory
    for (unsigned p = getD() + 1; p < getP(); ++p) {
        unsigned a = p;
        unsigned i_ = i(a);
        unsigned j_ = j(a);
        Sigma(i_,j_) = getParam(k, p);
#ifdef DETAIL_VERBOSE_2
        std::cerr << a << " " << i_ << " " << j_ << " " << k << " " << p <<" " << getParam(k,p) << std::endl;
#endif
    }
#ifdef DETAIL_VERBOSE_2
    std::cerr << Sigma;
#endif
    return Sigma;
}

double GaussianMixtureModelNDCommon::evalPDF(const unsigned k, const fvector_t& x) const {

    assert(m_theta . getThetas() . size() > 0);

    namespace ublas = boost::numeric::ublas;

    ublas::vector<double> zwe(getD()); // vector Zwischenergebnis
    ublas::vector<double> xmu = x - getMean(k);
    zwe = ublas::prod(getCachedInvSigma(k), xmu);

    double res = pow(2*M_PI, -0.5*(double)getD()) *
                    pow(getCachedSigmaDet(k), -0.5) * exp(- 0.5 *  ublas::inner_prod(xmu, zwe));

#ifdef DETAIL_VERBOSE_2
    if (isnan(res)) {
        std::cout << "Asked to evalPDF at sigma " << getSigmaMatrix(k) << ", mean " << getMean(k) << std::endl;
        std::cout << "(1) " << pow(2*M_PI, -0.5*(double)getD()) << " " << pow(getCachedSigmaDet(k), -0.5) << std::endl;
        std::cout << "(1) " << xmu << std::endl;
        std::cout << "(1) " << getCachedSigmaDet(k) << " " << getD() << std::endl;
        std::cout << "(1) " << zwe << std::endl;
        std::cout << "(2) " << exp(- 0.5 *  ublas::inner_prod(xmu, zwe)) << std::endl;
        std::cout << pow(2*M_PI, -0.5*(double)getD()) * pow(getCachedSigmaDet(k), -0.5) * exp(- 0.5 *  ublas::inner_prod(xmu, zwe)) << std::endl;
    }
#endif

    return res;
    // Seems OK
}

size_t GaussianMixtureModelNDCommon::getD() const {
#ifdef DETAIL_VERBOSE_2
        std::cerr << "GaussianMixtureModelNDCommon: getD() called.\n" << std::flush;
#endif
    return getDataObj() -> getDataDimensionality();
}
