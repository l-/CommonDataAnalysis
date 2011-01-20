/**
* @file GaussianMixtureModelNDParams.h++
* @version 0.141
* @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
*
*  Created on: Jan 20, 2011
*
*/

#pragma once

#include "EMTheta.h++"

#include <boost/numeric/ublas/symmetric.hpp>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

namespace CDA {

/**
* @class GaussianMixtureModelNDParams
*
* @brief A heteroscedastic N-dimensional Gaussian Mixture Model
*/
class GaussianMixtureModelNDParams : public EMThetas {

    /**
    * @brief This class should also know its K.
    * @todo avoid this duplpication
    */
    const int K;

    /**
    * @brief Data dimensionality, again (because it determines the dimensionality of matrices etc.)
    */
    const int D;

public:

    typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> sym_mtx_t;

    /**
    * @brief Fresh Constructor.
    *
    * @param[in] K_
    * @param[in] D_
    *
    * @section Parameter space dimensionality:
    * P = 1 + D + D(D+1)/2 (class + mean + covariances)!
    *
    * @section WARNING
    * Be careful with the parameter order!
    */
    GaussianMixtureModelNDParams(const unsigned K_, const unsigned D_);

    GaussianMixtureModelNDParams(const GaussianMixtureModelNDParams& other);

    /**
    * @brief Get param index - D of covariance parameter sigma(i,j) ;-)
    * i.e. does not take into account the other parameters, so the
    * offset still has to be added. clear?
    *
    * @section NOTA BENE
    * <b> \f$i <= j\f$ always! </b>
    * So the upper triangular part of the matrix is stored.
    */
    inline unsigned a(const unsigned i, const unsigned j) const {
        if (i>j) { return a(j,i); } // stupid me
        return i*getD() + j - (i*(i+1))/2;
    }

    /**
    * @brief Construct names of parameters, e.g. for output
    */
    virtual const std::string paramName(const unsigned p) const;


public: // why shouldn't they? parameters are set from the outside now

    /**
     * @brief Get i index of covariance parameter ;-) VERT
     * W/ OFFSET
     *
     * BROKEN!!!!!!!!!
     */
    inline unsigned i(const unsigned ain, bool remove_offset=true, const unsigned iter = 0) const {
        const unsigned D = getD();
        unsigned a = remove_offset ? ain - D - 1 : ain;
        assert(iter<D);
        if (a>=D-iter) { return i(a-D+iter, false, iter+1); }
        else {
            return a+iter;
        }
    }

    /**
     * @brief Get j index of covariance parameter ;-) HORZ
     * W/ OFFSET
     *
     * BROKEN!!!!!!!!
     */
    inline unsigned j(const unsigned ain, bool remove_offset=true, const unsigned iter = 0) const {
        const unsigned D = getD();
        unsigned a = remove_offset ? ain - D - 1 : ain;
        assert(iter<D);
        if (a>=D-iter) { return j(a-D+iter, false, iter+1); }
        else {
            return iter;
        }
    }

protected:

    std::vector<sym_mtx_t> m_cached_invsigmas;

    std::vector<double> m_cached_sigmadet;

    /**
    * @brief Structured datatype! Much cleaner and easier and less error-prone than before ...
    */
    typedef boost::tuple<double, fvector_t, sym_mtx_t> pdfparams_t;

    std::vector<pdfparams_t> thetas;

    /**
    * @brief Recalculate
    *
    * @return det(S)
    */
    double getSigmaDet(const unsigned k) const;

    /**
    * @brief Recalculate and Get inverse of Sigma of Gaussian no. k
    *
    * @return a fresh matrix
    */
    const sym_mtx_t getInvSigma(const unsigned k) const;

public:

    /**
    * @brief Update the cached functions of PDF parameters
    * To be called after each modification (i.e. after M-step).
    */
    void updateCached();

    /**
    * @brief Extract the covariance matrix of Gaussian no. k
    *
    * @param[in] k class no.
    *
    * @return a ublas matrix
    */
    const sym_mtx_t getSigmaMatrix(const unsigned k) const;

    /**
    * @brief
    *
    * @return det(S), if up to date
    */
    double getCachedSigmaDet(const unsigned k) const;

    /**
    * @brief Get inverse of Sigma of Gaussian no. k
    */
    const sym_mtx_t&
    getCachedInvSigma(const unsigned k) const;

    /**
    * @brief Get Mean vector of Gaussian no. k
    */
    const fvector_t& getMean(const unsigned k) const;
    // const boost::numeric::ublas::vector_range<const fvector_t> getMean(const unsigned k) const;

    /**
    * @brief get this D
    */
    unsigned getD() const;

    /**
    * @brief get this K
    */
    unsigned getK() const;

    /**
     * @brief To initialize, some setters
     * Covariance matrix of class k
     *
     * @section WARNING
     * don't forget to run updateCached() after initializing all classes
     *
     * @param k: class number
     */
    void setSigmaMatrix(const unsigned k);

    /**
     * @brief To initialize, some setters
     * Mean vector of class k
     *
     * @section WARNING
     * don't forget to run updateCached() after initializing all classes
     *
     * @param k: class number
     */
    void setMean(const unsigned k);

    /**
     * @brief To initialize, some setters
     * Class probability
     *
     * @section WARNING
     * don't forget to run updateCached() after initializing all classes
     *
     * @param k: class number
     */
    void setProb(const unsigned k);
};

} // namespace
