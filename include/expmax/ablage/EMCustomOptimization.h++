///**
// * @file EMCustomOptimization.h++
// *
// * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
// *
// *  Created on: Jan 5, 2011
// *
// */
//
//#pragma once
//
//#include "EMGenericMixtureModelCore.h++"
//
//namespace CDA {
//
///**
// * @class EMMbyGradientDescent
// *
// * @brief EM whose M-step is done by gradient descent
// */
//class EMMbyGradientDescent : public EMGenericMixtureModelCore {
//
//protected:
//
//    /**
//     * @brief Call superclass constructor only.
//     */
//    EMMbyGradientDescent(const unsigned D_, const unsigned K_, const unsigned P_)
//    : EMGenericMixtureModelCore(D_, K_, P_) {}
//
//    /**
//     * @brief The quasi-MLE step (the not-so-easy part)
//     *
//     * Implemented here.
//     */
//    void improveClusterModelParameters();
//
//    /**
//     * @brief if there are any values which are needed later
//     * @todo refactor!
//     */
//    virtual void preparations_for_evalPDFderivP() = 0;
//
//    virtual double evalPDFderivP(const unsigned k, const datapoint_t& x, const unsigned p) const = 0;
//
//    /**
//     * @brief Now if we are able to evaluate the PDF, even better: calculate a gradient,
//     * we can actually do something.
//     *
//     * @section NOTE
//     * This is used in the M-step to improve on all parameters.
//     *
//     * @param[in] p 0<=p<P
//     *
//     * @return The component of the gradient \f$\frac{\partial L}{\partial p}\f$ <b>(nicht wahr?)</b>
//     *
//     */
//    double improvementOnTheta(const unsigned p) const;
//};
//
//} //ns
