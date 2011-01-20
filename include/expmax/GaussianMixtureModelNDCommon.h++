/**
 * @file GaussianMixtureModelNDCommon.h++
 * @version 0.12
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 5, 2011
 *
 */

#pragma once

#include "EMData.h++"
#include "EMGenericMixtureModelCore.h++"
#include "GaussianMixtureModelNDParams.h++"
//#include "EMCustomOptimization.h++"

namespace CDA {

/**
 * @class GaussianMixtureModelNDCommon
 *
 * @brief A heteroscedastic N-dimensional Gaussian Mixture Model
 * minus the parameter optimization part, which is left purevirtual
 * thus a replacement for the EMData<fvector_t> class.
 *
 * Is-a EM.
 */
class GaussianMixtureModelNDCommon : public EMGenericMixtureModelCore<GaussianMixtureModelNDParams> {

public:

    typedef EMGenericMixtureModelCore<GaussianMixtureModelNDParams>::data_t data_t;
    typedef GaussianMixtureModelNDParams theta_t;
    typedef GaussianMixtureModelNDParams::sym_mtx_t sym_mtx_t;

protected:

    using EM<data_t, theta_t>::getDataObj;
    using EM<data_t, theta_t>::getThetaObj;

public:

    /**
     * @brief Constructor. Calls several superclass ones
     *
     * @param[in] K_
     * @param[in] D_
     * @param[in] data
     * @param[in] theta
     *
     * @section Parameter space dimensionality
     * P = D + D(D+1)/2 (mean + covariances)
     *
     * @section WARNING
     * Be careful with the parameter order!
     *
     */
    GaussianMixtureModelNDCommon(const unsigned K_, const unsigned D_, const data_t& data, const theta_t& theta);

    /**
     * @brief Evaluate model PDF of cluster k
     *
     * @param k class no.
     * @param x feature vector \f$\vec{x}\f$
     *
     * @return \f$p(x,k\vert\theta) = \frac{1}{(2\pi)^{\frac{D}{2}}\left|\Sigma_k\right|}e^{\frac{-1}{2}(\vec{x}-\vec{\mu_k})^T\Sigma_k^{-1}(\vec{x}-\vec{\mu_k})}\f$
     *
     */
    double evalPDF(const unsigned k, const datapoint_t& x) const;

    /**
     * @brief yet another ...
     */
    size_t getD() const;

    void setSigma(const unsigned k, const sym_mtx_t& sigma) { getThetaObj() -> getModifySigma(k) = sigma; }
    void setMean(const unsigned k, const fvector_t& param) { getThetaObj() -> getModifyMean(k) = param; }
    void setClassProb(const unsigned k, const double prob) { getThetaObj() -> setClassProb(k, prob); }

    /**
     * @todo this is only because the public getThetaObj() of EM:: isn't accessible somehow ...
     * I'm in need of a C++ expert here.
     * @param k
     * @param sigma
     */
    const sym_mtx_t& getSigma(const unsigned k) const { return getThetaObj() -> getSigmaMatrix(k); }
    const double getCachedSigmaDet(const unsigned k) const { return getThetaObj() -> getCachedSigmaDet(k); }
    const fvector_t& getMean(const unsigned k) const { return getThetaObj() -> getMean(k); }
    const double getClassProb(const unsigned k) const { return getThetaObj() -> getClassProb(k); }
};


} // namespace
