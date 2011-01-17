///**
// * @file FitMultivariateMulticlassByEM.h++
// *
// * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
// *
// *  Created on: Jan 5, 2011
// *
// */
//
//#pragma once
//
//#include "FitMulticlassByEM.h++"
//
//namespace CDA {
//
//
///**
// * @brief It's been demoted to a typedef.
// */
//class FitMultivariateMulticlassByEM : public FitMulticlassByEM<fvector_t> {
//
//public:
//
//    /**
//     * @brief Define the datapoint_t
//     */
//    typedef fvector_t datapoint_t;
//
//public:
//
//    /**
//     * @brief Constructor
//     *
//     * @param[in] K_ this way round, because of default parameter in superclass
//     * @param[in] D_
//     */
//    FitMultivariateMulticlassByEM(const unsigned K_, const unsigned D_)
//    : FitMulticlassByEM<fvector_t>(K_, D_) {}
//    // , EM<fvector_t>(D_) {}
//
//    /**
//     * @brief Convenience getter for D, same format as getN() and getP()
//     */
//    unsigned getD() const;
//};
//
//} // ns
