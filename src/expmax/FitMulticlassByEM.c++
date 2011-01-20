/**
* @file FitMulticlassByEM.c++
* @version 0.131
* @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
*
*  Created on: Jan 17, 2011
*
*/

#include "expmax/ProbabilisticClustering.h++"

#include <boost/foreach.hpp>

#include <boost/mem_fn.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include "utils/comparisons.h++"

// Loop-avoiding trickery
#include <boost/iterator/counting_iterator.hpp>

using namespace CDA;

template<class data_T, class theta_T>
const std::string FitMulticlassByEM<data_T, theta_T>::dumpParameters(const unsigned k, const boost::optional<unsigned int> iteration) const {
    std::stringstream out;

    if (iteration)
        out << *iteration << ";";

    out << k << ";";

    unsigned i=0;
    for (unsigned n=0; n<getThetaObj() -> getP(); ++n) {
        double vi = getThetaObj() -> getThetas(k,n);
        out << vi;

        if (i < getK() - 1) {
            out << ";";
        } else {
            out << "\n";
        }

        ++i;
    }

    return out.str();
}

template<class data_T, class theta_T>
const std::string FitMulticlassByEM<data_T, theta_T>::dumpParameters(const boost::optional<unsigned int> iteration) const {
    std::stringstream out;
    for (unsigned i=0; i<K; ++i)
        out << dumpParameters(i, iteration);
    return out.str();
}

template<class data_T, class theta_T>
double FitMulticlassByEM<data_T, theta_T>::logLikelihood() const {
    double res = 0;
    for (unsigned n=0; n<getN(); ++n) {
        double beitraege = 0;
        for (unsigned k=0; k<getK(); ++k) {
            beitraege += getPk(k) * evalPDF(k, getDataObj()->getData(n));
        }
        res += log(beitraege);
    }
    return res;
}

template<class data_T, class theta_T>
const fvector_t&  FitMulticlassByEM<data_T, theta_T>::getHiddenParamEstimate(const unsigned n) const {
    return classif[n];
}

template<class data_T, class theta_T>
unsigned int FitMulticlassByEM<data_T, theta_T>::getBestClass(const unsigned n) const {
    const fvector_t& cs = getHiddenParamEstimate(n);
    return *( std::max_element<boost::counting_iterator<int>, compareByVectorElementValue<fvector_t, size_t> >
    (boost::counting_iterator<int>(0),
            boost::counting_iterator<int>(getK()),
            compareByVectorElementValue<fvector_t>(cs)) );
}

template<class data_T, class theta_T>
std::vector<fvector_t>& FitMulticlassByEM<data_T, theta_T>::getModifyClassif() {
    return classif;
}

template<class data_T, class theta_T>
const std::vector<fvector_t>& FitMulticlassByEM<data_T, theta_T>::getClassif() const {
    return classif;
}

template<class data_T, class theta_T>
const fvector_t& FitMulticlassByEM<data_T, theta_T>::getClassif(const unsigned n) const {
    return classif[n];
}

template<class data_T, class theta_T>
std::vector<unsigned int> FitMulticlassByEM<data_T, theta_T>::getClassifList() const {
    std::vector<unsigned int> result;
    std::for_each(boost::counting_iterator<int>(0),
            boost::counting_iterator<int>(getN()),
            boost::lambda::bind(boost::mem_fn(&std::vector<unsigned int>::push_back),
                    &result,
                    boost::lambda::bind(
                            boost::mem_fn(&FitMulticlassByEM<data_T, theta_T>::getBestClass),
                            this,
                            boost::lambda::_1)));
    return result;
}

template<class data_T, class theta_T>
void FitMulticlassByEM<data_T, theta_T>::initClassif() {
    classif . clear();

#ifdef VERBOSE
    std::cerr << "initClassif called\n";
#endif

    for (unsigned i=0; i<getN(); ++i) {
        classif . push_back ( fvector_t(K, 1/(double)K) );
    }
}


template<class data_T, class theta_T>
double FitMulticlassByEM<data_T, theta_T>::getPk(const unsigned k) const {
    return getThetaObj() -> getThetas(k, 0);
}

template<class data_T, class theta_T>
void FitMulticlassByEM<data_T, theta_T>::update_hidden() {

#ifdef VERBOSE
    std::cout << className() << ": E step.\n";
#endif

    // @todo triple-check everything

    for (unsigned n=0; n<getN(); ++n) {
        for (unsigned k=0; k<getK(); ++k) { // m_theta.getThetas().size()

#ifdef DETAIL_VERBOSE_2
            // maybe only the first time
             std::cout << "Paramset number " << k << " -- " << getPk(k) << " " << evalPDF(k, getDataObj()->getData(n)) << " " << classif[n](k) << std::endl;
#endif
            classif[n](k) = getPk(k) * evalPDF(k, getDataObj()->getData(n));
            if (isnan(classif[n](k))) {
                classif[n](k) = 0.0;
#ifdef VERBOSE_2
                // maybe only the first time
            std::cerr << "Error: encountered NaN\n";
#endif
            } // wtf?
            // @todo better error handling
        }

        if (norm_1(classif[n]) != 0)
            classif[n] /= norm_1(classif[n]);

        //                            classif[n] = FLT_EPSILON; // make sense?
    }
}

template<class data_T, class theta_T>
double FitMulticlassByEM<data_T, theta_T>::getParam(const unsigned k, const unsigned p) const {
    return getThetaObj() -> getThetas(k,p);
}

template<class data_T, class theta_T>
unsigned int FitMulticlassByEM<data_T, theta_T>::getK() const
{
    return K;
}

// Without these two declarations, the file would be useless ;-)
template class FitMulticlassByEM<EMData<double>, EMThetas>;
template class FitMulticlassByEM<VectorEMData, EMThetas>;
template class FitMulticlassByEM<VectorEMData, GaussianMixtureModelNDParams>;
