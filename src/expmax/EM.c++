/**
 * @file EM.c++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 */

#include "expmax/EM.h++"

using namespace CDA;

// => EMData itself.
//
//// templates may not be virtual ...
//template<class datapoint_t>
//EMData<datapoint_t>& EM<datapoint_t>::getDataObj() {
//    return m_data;
//}
//
//template<class datapoint_t>
//const EMData<datapoint_t>& EM<datapoint_t>::getDataObj() const {
//    return m_data;
//}

template<class datapoint_t>
unsigned int EM<datapoint_t>::getN() const {
    return getDataObj() . getNumberOfDataPoints();
}

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

template<class datapoint_t>
EM<datapoint_t>::EM()
  {  }

template class EM<double>;
template class EM<fvector_t>;
