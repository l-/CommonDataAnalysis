/**
 * @file ProbabilisticClustering.c++
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 * @version 0.1
 *
 * @section LICENSE
 *
 * This file is released under LGPL v3.0.
 *
 *  Created on: Dec 16, 2010
 *
 *  @todo crippling bugs in multivariate gaussian specialization
 *
 */

#include "ProbabilisticClustering.h++"

#include <sstream>
#include <boost/numeric/ublas/triangular.hpp>

using namespace CDA;

// template<>
// size_t EMData<double>::getDataDimensionality() const { return 1; }
// template<>
// size_t EMData<fvector_t>::getDataDimensionality() const { return D; }

namespace CDA {

template<> EMData<double>::EMData(const unsigned D_) : D(1) {}; // all's well
template<> EMData<fvector_t>::EMData(const unsigned D_) : D(D_) {};

} // ns

template class EMData<double>;
template class EMData<fvector_t>;

template<class datapoint_t>
const std::string FitMulticlassByEM<datapoint_t>::dumpParameters(const unsigned k, const boost::optional<unsigned int> iteration) const {
    std::stringstream out;

    if (iteration)
        out << *iteration << ";";

    out << k << ";";

    unsigned i=0;
    BOOST_FOREACH ( double vi, m_theta.getThetas()[k] ) {
        out << vi;

        if (i < m_theta.getThetas()[k].size() - 1) {
            out << ";";
        } else {
            out << "\n";
        }

        ++i;
    }

    return out.str();
}

template<class datapoint_t>
const std::string FitMulticlassByEM<datapoint_t>::dumpParameters(const boost::optional<unsigned int> iteration) const {
    std::stringstream out;
    for (unsigned i=0; i<K; ++i)
        out << dumpParameters(i, iteration);
    return out.str();
}

template<class datapoint_t>
double FitMulticlassByEM<datapoint_t>::logLikelihood() const {
    double res = 0;
    for (unsigned n=0; n<m_data.getN(); ++n) {
        double beitraege = 0;
        for (unsigned k=0; k<m_theta.getThetas().size(); ++k) {
            beitraege += getPk(k) * evalPDF(k, m_data.getData(n));
        }
        res += log(beitraege);
    }
    return res;
}


unsigned FitMultivariateMulticlassByEM::getD() const {
    return m_data . getDataDimensionality();
}

template<class datapoint_t>
double FitMulticlassByEM<datapoint_t>::getPk(const unsigned k) const {
    return m_theta.getThetas()[k](0);
}

template<class datatype_t>
void FitMulticlassByEM<datatype_t>::update_hidden() {

#ifdef VERBOSE
    std::cout << className() << ": E step.\n";
#endif

    // @todo triple-check everything

    for (unsigned n=0; n<m_data.getN(); ++n) {
        for (unsigned k=0; k<K; ++k) { // m_theta.getThetas().size()
            // std::cout << "Paramset number " << k << " -- " << getPk(k) << " " << evalPDF(k, data[n]) << " " << classif[n](k) << std::endl;
            classif[n](k) = getPk(k) * evalPDF(k, m_data.getData(n));
            if (isnan(classif[n](k))) { classif[n](k) = 0.0;
#ifdef EXTRA_VERBOSE // maybe only the first time
                std::cerr << "Error: encountered NaN\n";
#endif
            } // wtf?
            // @todo better error handling
        }

        classif[n] /= norm_1(classif[n]);

        //                            classif[n] = FLT_EPSILON; // make sense?
        //                        }

        // std::cout << data[n] << "as" << classif[n] << std::endl;
    }
}

template<class datapoint_t>
double FitMulticlassByEM<datapoint_t>::getParam(const unsigned k, const unsigned p) const {
    return m_theta.getThetas()[k][p+1];
}

template class FitMulticlassByEM<double>;
template class FitMulticlassByEM<fvector_t>;

template<class datapoint_t>
void EM<datapoint_t>::EMrun(const unsigned MAXITER, const double thresh, const boost::optional<std::ostream*> output_csv) {

    assert(m_theta . getThetas() . size() > 0);

    // @todo generic output-handler

    double q_likelihood = -INFINITY, q_likelihood_old = -INFINITY;
    double improvement = 0;

    if (output_csv)
        (**output_csv) << getCSVHeader();

    if (output_csv) { **output_csv << dumpParameters(boost::make_optional((unsigned int)0)) << std::flush; }

    // How to select the threshold??? Likelihoods can be very small ...

    for (unsigned c=1; c<MAXITER; ++c) {

        q_likelihood_old = q_likelihood;
        EMstep();
        q_likelihood = logLikelihood();

        improvement = q_likelihood - q_likelihood_old;

        std::cout << q_likelihood << std::endl;

        if (output_csv)
            { **output_csv << dumpParameters(boost::make_optional(c)) << std::flush; }

        if (improvement < thresh) {
            break;
        }
    }

    std::cout << "Done, with log-likelihood " << q_likelihood << std::endl << std::flush;
}

void GaussianMixtureModel1D::update_thetas() {

#ifdef VERBOSE
    std::cout << className() << ": M step.\n";
#endif

    double sumclassif[K];
    for (unsigned k=0; k<K; ++k) {
        sumclassif[k] = 0.0;
    }

    // Use thetas[k] as accumulators, since in effect the estimates are based on simple sample statistics.
    for (unsigned k=0; k<K; ++k) {
        m_theta.getModifyThetas()[k](2) = 0;
        m_theta.getModifyThetas()[k](1) = 0;
        m_theta.getModifyThetas()[k](0) = 0;
    }

    for (unsigned n=0; n<m_data.getN(); ++n) {
        for (unsigned k=0; k<K; ++k) {
            sumclassif[k] += classif[n](k);
            m_theta.getModifyThetas()[k](1) += (classif[n](k) * m_data.getData(n));
        }
    }

    for (unsigned k=0; k<K; ++k) {
        m_theta.getModifyThetas()[k](1) /= sumclassif[k];
    }

    for (unsigned n=0; n<m_data.getN(); ++n) {
        for (unsigned k=0; k<K; ++k) {
            // using NEW mean, thus in extra step.

            // @todo: effizienter und numerisch besser berechnen ... aber erstmal das ganze verfahren geradebiegen.
            m_theta.getModifyThetas()[k](2) +=
                classif[n](k) * squaredDistanceToMean(k, m_data.getData(n));
        }
    }

    for (unsigned k=0; k<K; ++k) {
        m_theta.getModifyThetas()[k](2) = sqrt(m_theta.getModifyThetas()[k](2) / sumclassif[k]);
        m_theta.getModifyThetas()[k](0) = 1/((double)(m_data.getN())) * sumclassif[k];
    }

    // @todo atomic
}

const std::string GaussianMixtureModel1D::className() const {
    return std::string("GaussianMixtureModel1D");
}

const std::string GaussianMixtureModel1D::getCSVHeader() const {
    return std::string("iter;k;p;m;sigma\n"); // @todo variable num classes
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

const std::string CircularMixtureModel1D::className() const {
    return std::string("CircularMixtureModel1D");
}

const std::string CircularMixtureModel1D::getCSVHeader() const {
    return std::string("iter;k;p;mu;kappa\n"); // @todo variable num classes
}

double CircularMixtureModel1D::getMu(const unsigned k) const {
    return m_theta.getThetas()[k](1);
}

double CircularMixtureModel1D::getKappa(const unsigned k) const {
    return m_theta.getThetas()[k](2);
}

double CircularMixtureModel1D::evalPDF(const unsigned k, const double& x) const {
    return VonMises()(getMu(k), getKappa(k), x);
}

void CircularMixtureModel1D::update_thetas() {

    // @todo confirm this is indeed correct

#ifdef VERBOSE
    std::cout << className() << ": M step.\n";
#endif

    double sumclassif[K];
    for (unsigned k=0; k<K; ++k) {
        sumclassif[k] = 0.0;
    }

    // Use thetas[k] as accumulators, since in effect the estimates are based on simple sample statistics.
    for (unsigned k=0; k<K; ++k) {
        m_theta.getModifyThetas()[k](2) = 0;
        m_theta.getModifyThetas()[k](1) = 0;
        m_theta.getModifyThetas()[k](0) = 0;
    }

    // Finish estimate of class probability
    for (unsigned k=0; k<K; ++k) {
        m_theta.getModifyThetas()[k](0) = 1/((double)(m_data.getN())) * sumclassif[k];
    }

    // Estimate of mu

    // @todo FINISH THIS!!! DEC28
}

const std::string EMGenericMixtureModelCore::className() const {
    return std::string("EMGenericMixtureModelCore"); // shouldn't occur, class is still abstract
}

void EMGenericMixtureModelCore::update_thetas() {

#ifdef VERBOSE
    std::cout << className() << ": M step.\n";
#endif

    // The common part: updating of overall class probabilities
    double sumclassif[K];
    for (unsigned k=0; k<K; ++k) {
        sumclassif[k] = 0.0; }

    for (unsigned n=0; n<m_data.getN(); ++n) {
        for (unsigned k=0; k<K; ++k) {
            sumclassif[k] += classif[n](k);
        }
    }

    for (unsigned k=0; k<K; ++k) {
        m_theta.getModifyThetas()[k](0) = 1/((double)(m_data.getN())) * sumclassif[k]; }

    improveClusterModelParameters();
}

const std::string EMGenericMixtureModelCore::getCSVHeader() const {
    std::stringstream out;

    out << "iteration" << ";";
    out << "k" << ";";
    out << "p" << ";";
    for (unsigned p=0; p<P; ++p) {
        out << paramName(p);
        if (p < P - 1) {
            out << ";";
        } else {
            out << "\n";
        }
    }

    return out.str();
}

void EMMbyGradientDescent::improveClusterModelParameters() {

    // Timidly trying to move upwards. Do a line search? Implement a real optimization method?
    // @todo The derivatives have to go. Do it numerically.
    const double GAMMA = 0.1;

    for (int c=0; c<50; ++c) {

        preparations_for_evalPDFderivP();

        for (unsigned k = 0; k < getK(); ++k) {

            // evalPDF(k)
            fvector_t improvement(getP() + 1);

            improvement(0) = 0; // hackish, @todo further beautify internal interfaces

            for (unsigned p = 0; p < getP(); ++p) {
                improvement(p + 1) = improvementOnTheta(p);
            }

#ifdef VERBOSE
            std::cerr << "\"IMPROVEMENT\"   " << improvement << std::endl;
            std::cerr << "\"LogLikelihood\" " << logLikelihood() << std::endl;
            std::cerr << dumpParameters(k) << std::endl;
#endif


            if (boost::numeric::ublas::norm_2(improvement) < 1e-5) break;

            m_theta . getModifyThetas()[k] += GAMMA * improvement / boost::numeric::ublas::norm_2(improvement);
        }
    }
}

double EMMbyGradientDescent::improvementOnTheta(const unsigned p) const {

    double result = 0.0;

    // For all feature points. like in the expression of the likelihood, which we want to maximize here.

    for (unsigned n=0; n<m_data.getN(); ++n) {

        double sum_prob = 0.0;
        double sum_der = 0.0;

        // Is it really that easy? The formulas say so. We'll see ...

        for (unsigned k=0; k<K; ++k) {

            sum_prob += getPk(k) * evalPDF(k, m_data.getData(n));
            sum_der  += getPk(k) * evalPDFderivP(k, m_data.getData(n), p);

        // if (n==0)
//                std::cerr << k << " " << p << " " << evalPDF(k, m_data.getData(n)) << " " << evalPDFderivP(k, m_data.getData(n), p) << std::endl;
        }

//        std::cerr << sum_prob << " " << sum_der << std::endl;

        result += (sum_der / sum_prob);
    }

    return (isnan(result) || isinf(result) ? 0.0 : result / (double)m_data.getN()); // Yes???
}

const std::string GaussianMixtureModel::paramName(const unsigned p) const {
    std::stringstream out;
    if (p <= m_data . getDataDimensionality()) {
        out << "m_" << p;
    } else {
        int a = p - m_data . getDataDimensionality();
        out << "s_" << i(a) << "_" << j(a) << std::endl;
    }
    return out.str();
}

const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> GaussianMixtureModel::getSigmaMatrix(const unsigned k) const {
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

double GaussianMixtureModel::getSigmaDet(const unsigned k) const {
    return det(getSigmaMatrix(k));
}

const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper>
GaussianMixtureModel::getInvSigma(const unsigned k) const {
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

const boost::numeric::ublas::vector_range<const fvector_t>
GaussianMixtureModel::getMean(const unsigned k) const {
    using namespace boost::numeric::ublas;
    return project(m_theta.getThetas()[k], range(1, 1+getD()));
}

double GaussianMixtureModel::evalPDF(const unsigned k, const datapoint_t& x) const {

    assert(m_theta . getThetas() . size() > 0);

    namespace ublas = boost::numeric::ublas;

    ublas::vector<double> zwe(getD()); // Zwischenergebnis
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

double GaussianMixtureModel::getSigmaCofactor(const unsigned k, const unsigned i, const unsigned j) const {
    // @todo make it better
    boost::numeric::ublas::matrix<double> mtx(getD()-1, getD()-1);
    boost::numeric::ublas::matrix<double> sigma(getSigmaMatrix(k));
    // Now this is completely crappy. The whole module has degenerated to a mere proof of concept.

    for (unsigned j_=0; j_<getD(); ++j_) {
        for (unsigned i_=0; i_<getD(); ++i_) {

            if (i_ != i && j_ != j) {
                mtx(i_>=i ? i_-1 : i_,
                    j_>=j ? j_-1 : j_) = sigma(i_, j_);
            }
        }
    }

    // std::cerr<< "Cofactor " << i << " " << j << "  of " << sigma << " is " << mtx <<    std::endl;

    return ( (i+j)%2==0 ? 1 : -1 ) * det(mtx);
}

void GaussianMixtureModel::preparations_for_evalPDFderivP() {
    m_cached_invsigmas.clear();
    for (unsigned k=0; k<getK(); ++k) {
        m_cached_invsigmas . push_back( getInvSigma(k) );
    }
}

double GaussianMixtureModel::getInvSigmaDerivTerm(const unsigned k, const unsigned i, const unsigned j, const datapoint_t& x) const {
    namespace ublas = boost::numeric::ublas;
    ublas::vector<double> zwe(getD());
    zwe = ublas::prod(getInvSigma(k), x - getMean(k));

    return zwe(i) * zwe(j); // now that's more like it.
}

double GaussianMixtureModel::evalPDFderivP(const unsigned k, const datapoint_t& x, const unsigned p) const {

    namespace ublas = boost::numeric::ublas;

    if (p < getD()) {
        // easy
        ublas::vector<double> transxm(getD());
        transxm = ublas::prod(getInvSigma(k), x - getMean(k));
        // @todo: please quintuple-check every expression

        return 2 * transxm(p) * evalPDF(k, x);

    } else {
        unsigned a = p - getD();
        unsigned i_ = i(a);
        unsigned j_ = j(a);

        double result = evalPDF(k, x) * ( -1/(2*getSigmaDet(k)) * getSigmaCofactor(k, i_, j_) - 1/2 * getInvSigmaDerivTerm(k, i_, j_, x) );

        return result;
    }
}

template class EM<double>;
template class EM<fvector_t>;
