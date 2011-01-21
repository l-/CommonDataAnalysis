/**
 * @file CircularStatistics.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Dec 28, 2010
 *
 * @brief
 *  Special functions for circular statistics
 */

#pragma once

// Special functions from GSL
#include <gsl/gsl_sf_bessel.h>
#include <cmath>

// Hell, just fake it
double wrap2pi(const double phi);

namespace CDA {

/**
 * Von Mises PDF and derivatives for circular statistics
 *
 * (it approximates a wrapped Gaussian for small \f$Kappa\f$
 *  and becomes a uniform distribution for large \f$Kappa\f$).
 */
struct VonMises {

    /**
     * Get the function value
     *
     * \f[
     *   M(\mu, \kappa, x) = \frac{\exp{\kappa \cdot cos (x-\mu)}}{2\pi I_0 (\kappa)}
     * \f]
     */
    inline double operator()(const double mu, const double kappa, const double phi) const {
        return exp(kappa * cos(phi - mu)) / (2 * M_PI * gsl_sf_bessel_I0 ( kappa ));
    }

    // @todo should be static. 

    /**
     * Get the value of \f$\frac{\partial M}{\partial \mu}\f$
     *
     * Important for parameter estimation
     */
    inline
    double deriv_mu(const double mu, const double kappa, const double phi) const {
        return kappa * exp(kappa * cos(phi - mu)) * sin(phi - mu) / (2 * M_PI * gsl_sf_bessel_I0 ( kappa ));
    }

    /**
     * Get the value of \f$\frac{\partial M}{\partial \kappa}\f$
     *
     * Important for parameter estimation
     */
    inline
    double deriv_kappa(const double mu, const double kappa, const double phi) const {
        const double f  = exp(kappa * cos(phi - mu));
        const double f_ = cos(phi - mu) * exp(kappa * cos(phi - mu));
        const double g  = 2 * M_PI * gsl_sf_bessel_I0 ( kappa );
        const double g_ = 2 * M_PI * gsl_sf_bessel_In ( 1, kappa ); // true?
        return (f_*g-f*g_)/(g*g);

        // = M * (-I1/I0 + cos(phi - mu)) ?
    }

//    double m_mu;
//    double m_kappa;
//
//    VonMises(const double mu, const double kappa)
//      : m_mu(mu)
//      , m_kappa(kappa) { }
//
//    inline double operator()(const double phi) const {
//        return operator()(m_mu, m_kappa, phi);
//    }

};


} // namespace
