/**
* @file FitMulticlassByEM.c++
*
* @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
*
*  Created on: Jan 17, 2011
*
*/

#include "expmax/FitMulticlassByEM.h++"

#include <boost/foreach.hpp>

#include <boost/mem_fn.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include "utils/comparisons.h++"

#include <boost/iterator/counting_iterator.hpp>


using namespace CDA;

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
    for (unsigned n=0; n<getN(); ++n) {
        double beitraege = 0;
        for (unsigned k=0; k<m_theta.getThetas().size(); ++k) {
            beitraege += getPk(k) * evalPDF(k, getDataObj().getData(n));
        }
        res += log(beitraege);
    }
    return res;
}

template<class datapoint_t>
const fvector_t&  FitMulticlassByEM<datapoint_t>::getHiddenParamEstimate(const unsigned n) const {
    return classif[n];
}

template<class datapoint_t>
unsigned int FitMulticlassByEM<datapoint_t>::getBestClass(const unsigned n) const {
    const fvector_t& cs = getHiddenParamEstimate(n);
    return *( std::max_element<boost::counting_iterator<int>, compareByVectorElementValue<fvector_t, size_t> >
    (boost::counting_iterator<int>(0),
            boost::counting_iterator<int>(getK()),
            compareByVectorElementValue<fvector_t>(cs)) );
}

template<class datapoint_t>
std::vector<unsigned int> FitMulticlassByEM<datapoint_t>::getClassifList() const {
    std::vector<unsigned int> result;
    std::for_each(boost::counting_iterator<int>(0),
            boost::counting_iterator<int>(getN()),
            boost::lambda::bind(boost::mem_fn(&std::vector<unsigned int>::push_back),
                    &result,
                    boost::lambda::bind(
                            boost::mem_fn(&FitMulticlassByEM<datapoint_t>::getBestClass),
                            this,
                            boost::lambda::_1)));
    return result;
}

template<class datapoint_t>
void FitMulticlassByEM<datapoint_t>::initClassif() {
    classif . clear();

    for (unsigned i=0; i<getN(); ++i) {
        classif . push_back ( fvector_t(K, 1/(double)K) );
    }
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

    for (unsigned n=0; n<getN(); ++n) {
        for (unsigned k=0; k<getK(); ++k) { // m_theta.getThetas().size()

#ifdef EXTRA_VERBOSE
            // maybe only the first time
             std::cout << "Paramset number " << k << " -- " << getPk(k) << " " << evalPDF(k, getDataObj().getData(n)) << " " << classif[n](k) << std::endl;
#endif
            classif[n](k) = getPk(k) * evalPDF(k, getDataObj().getData(n));
            if (isnan(classif[n](k))) {
                classif[n](k) = 0.0;
#ifdef VERBOSE
                // maybe only the first time
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

template<class datapoint_t>
unsigned int FitMulticlassByEM<datapoint_t>::getK() const
{
    return K;
}


// Without these two declarations, the file would be useless ;-)
template class FitMulticlassByEM<double>;
template class FitMulticlassByEM<fvector_t>;

