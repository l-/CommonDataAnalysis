/**
 * @file EM.c++
 * @version 0.12
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 */

#include "expmax/EM.h++"

// Common instantiations, if not in this file the code will not be compiled
#include "expmax/VectorEMData.h++"
#include "expmax/GaussianMixtureModelNDParams.h++"

using namespace CDA;

template<class data_t, class theta_t>
unsigned int EM<data_t, theta_t>::getN() const {
    return getDataObj() -> getNumberOfDataPoints();
}

template<class data_t, class theta_t>
void EM<data_t, theta_t>::EMrun(const unsigned MAXITER, const double thresh, const boost::optional<std::ostream*> output_csv) {

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

template<class data_t, class theta_t>
EM<data_t, theta_t>::EM(const data_t& data, const theta_t& theta) : m_data(data), m_theta(theta)
  {  }
template<class data_t, class theta_t>
EM<data_t, theta_t>::EM(const EM& other)
  : m_data(other.m_data)
  , m_theta(other.m_theta) {}

template<class data_t, class theta_t>
data_t* EM<data_t, theta_t>::getDataObj() { return &m_data; }
template<class data_t, class theta_t>
const data_t* EM<data_t, theta_t>::getDataObj() const { return &m_data; }
template<class data_t, class theta_t>
theta_t* EM<data_t, theta_t>::getThetaObj() { return &m_theta; }
template<class data_t, class theta_t>
const theta_t* EM<data_t, theta_t>::getThetaObj() const { return &m_theta; }

template class EM<EMData<double>, EMThetas>;
template class EM<VectorEMData, EMThetas>;
template class EM<VectorEMData, GaussianMixtureModelNDParams>;
