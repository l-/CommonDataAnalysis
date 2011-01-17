/**
 * @file EMGenericMixtureModelCore.h++
 *
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
 *
 * @section NOTE
 * for multivariate.
 */
class EMGenericMixtureModelCore : public FitMultivariateMulticlassByEM {

public:
    /**
     * @brief Propagate the datapoint_t
     */
    typedef FitMultivariateMulticlassByEM::datapoint_t datapoint_t;

protected:

    using FitMultivariateMulticlassByEM::classif;
    using FitMultivariateMulticlassByEM::getN;

    /**
     * @brief Parameter space dimensionality.
     *
     * @section NOTA BENE
     * We don't count the class probability here,
     * as it is handled differently.
     */
    const unsigned P;

    /**
     * @brief <b>M-step</b>: update parameters of model PDF
     *
     * It calls improveClusterModelParameters, where the implementer
     * specifies how the other parameters (beside a priori class probabilities)
     * are updated.
     */
    void update_thetas();

    /**
     * @brief The quasi-MLE step. Since this may be done in different-different ways,
     * it can be mixed in from a separate class.
     *
     * It is called by update_thetas
     */
    virtual void improveClusterModelParameters() = 0;

    /**
     * @brief Get param space dimensionality.
     *
     * @return If it is still zero, there's a programming error!
     */
    unsigned int getP() const { assert(P!=0); return P; }

public:

    /**
     * @brief Constructor to be called by implementer, please
     *
     * @param[in] K_
     * @param[in] P_ specific to the model.
     */
    EMGenericMixtureModelCore(const unsigned K_, const unsigned P_)
      : FitMultivariateMulticlassByEM(K_), P(P_) {}

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
     * @param[in] p 0<=p<P
     */
    virtual const std::string paramName(const unsigned p) const = 0;

    /**
     * @brief Each class gets its own line, that's the easiest way
     * Assume iteration no. is to be shown always
     *
     * @return String the parameter names together.
     */
    const std::string getCSVHeader() const;

};


} //ns
