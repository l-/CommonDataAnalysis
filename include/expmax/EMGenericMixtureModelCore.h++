/**
 * @file EMGenericMixtureModelCore.h++
 * @version 0.12
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 5, 2011
 *
 */

#pragma once

#include "FitMultivariateMulticlassByEM.h++"

namespace CDA {

/**
 * @class EMGenericMixtureModelCore
 *
 * @brief Ideally, each model knows its data layout, parameter estimators and stuff.
 * There could be a general blueprint for dealing with models which cannot be handled
 * analytically, prompting a gradient search for the MLE step
 */
template <class theta_T>
class MYEXPORT EMGenericMixtureModelCore : public FitMultivariateMulticlassByEM<theta_T> {

public:
    /**
     * @brief Propagate the datapoint_t
     */
    typedef typename FitMultivariateMulticlassByEM<theta_T>::data_t data_t;
    typedef fvector_t datapoint_t;
    typedef theta_T theta_t;

    /**
     * @brief C++ gurus: Why should this be needed?
     */
    using FitMulticlassByEM<data_t, theta_t>::getK;
    using EM<data_t, theta_t>::getN;

    /**
     * @brief Which is protected ...
     */
    using EM<data_t, theta_t>::getDataObj; // which is virtual, rather that m_data

    /**
     * @brief same
     */
    using EM<data_t, theta_t>::getThetaObj;

    /**
     * @brief similar
     */
    using FitMultivariateMulticlassByEM<theta_t>::getClassif;

protected:

    using FitMultivariateMulticlassByEM<theta_T>::getModifyClassif;

    /**
     * @brief Parameter space dimensionality.
     *
     * @section NOTA BENE
     * We didn't count the class probability before,
     * as it was handled differently.
     * ONCE AND FOR ALL (but not for all times),
     * P is counted INCLUDING the class prior probability.
     */
    const unsigned P;

    /**
     * @brief <b>M-step</b>: update parameters of model PDF
     *
     * @section DESCRIPTION
     * It calls improveClusterModelParameters, where the implementer
     * specifies how the other parameters (beside a priori class probabilities)
     * are updated.
     */
    void update_thetas();

    /**
     * @brief The quasi-MLE step. Since this may be done in different-different ways,
     * it can be mixed in from a separate class.
     *
     * @section DESCRIPTION
     * It is called by update_thetas
     */
    virtual void improveClusterModelParameters() = 0;

public:

    /**
     * @brief Constructor to be called by implementer, please
     *
     * @param[in] K_
     * @param[in] P_ specific to the model.
     * @param[in] D_ // really belongs here? don't think so.
     * @param[in] data
     * @param[in] theta
     */
    EMGenericMixtureModelCore(const unsigned K_, const unsigned P_, const unsigned D_, const data_t& data, const theta_t& theta)
      : FitMultivariateMulticlassByEM<theta_t>(K_, data, theta), P(P_) {
#ifdef VERBOSE_2
        std::cerr << "EMGenericMixtureModelCore Constructor called, P should equal " << P << std::endl;
#endif
        }

    /**
     * @brief The individual way of evaluating the PDF
     *
     * @param k class no.
     * @param x
     *
     * @return probability \f$p(x,k\vert\theta)\f$
     */
    virtual double evalPDF(const unsigned k, const datapoint_t& x) const = 0;

    /**
     * @brief Class name for logging ...
     */
    virtual const std::string className() const;

    /**
     * @brief Name of a parameter p=\f$\theta_p\f$
     * (For output)
     *
     * avoid this variant; redefine it.
     *
     * @param[in] p 0<=p<P
     */
    virtual const std::string paramName(const unsigned p) const {
        std::stringstream sstr; sstr<<"p"<<p;
        return sstr.str();
    }

    /**
     * @brief Each class gets its own line, that's the easiest way
     * Assume iteration no. is to be shown always
     *
     * @return String the parameter names together.
     */
    const std::string getCSVHeader() const;

    /**
     * @brief Get param space dimensionality.
     *
     * @return If it is still zero, there's a programming error!
     */
    unsigned int getP() const { assert(P!=0); return P; }

};

} //ns
