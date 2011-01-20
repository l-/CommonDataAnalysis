/**
 * @file GaussianMixtureModelNDParams.c++
 * @version 0.141
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 20, 2011
 *
 */

#include "GaussianMixtureModelNDParams.h++"

using namespace CDA;

const GaussianMixtureModelNDParams::sym_mtx_t
GaussianMixtureModelNDParams::getSigmaMatrix(const unsigned k) const {
    return thetas[k].get<2>();
}

unsigned GaussianMixtureModelNDParams::getD() const { return D; }
unsigned GaussianMixtureModelNDParams::getK() const { return K; }

GaussianMixtureModelNDParams::GaussianMixtureModelNDParams(const unsigned K_, const unsigned D_)
  : K(K_), D(D_) {
    // It shall never be in an undefined state!
  for (int k=0; k<getK(); ++k) {
      thetas.push_back(pdfparams_t(1/(double)k,
                       boost::numeric::ublas::zero_vector<double>(D),
                       boost::numeric::ublas::identity_matrix<double>(D)));
  }
}

GaussianMixtureModelNDParams::GaussianMixtureModelNDParams(const GaussianMixtureModelNDParams& other);
