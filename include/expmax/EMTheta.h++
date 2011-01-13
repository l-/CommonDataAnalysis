/**
 * @file EMTheta.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 5, 2011
 *
 */

#pragma once
#include "common_definitions.h++"

namespace CDA {

/**
 * @class EMThetas
 *
 * @brief Model parameters for EM
 *
 * @section DESCRIPTION
 * To avoid duplication. Can't just put it into EM and inherit
 * because of duplicate class problems
 *
 * @section ARCHITECTURE CAVEAT
 * I dislike having this class around, this should all be in EM
 * in any proper programming language, it would be possible.
 */
class EMThetas {

protected:
    /**
     * @brief Model parameters
     * possibly a single vector, or K PDF parameter vectors, one per class
     * => these belong in the concrete classes? I hate duplication!
     */
    std::vector<fvector_t> thetas;

public:
    /**
     * Initial guess for model parameters -- strictly needed
     *
     * Input data format differs for each individual distribution,
     * be careful and read the documentation.
     *
     * @todo random numbers when none given
     */
    template<class II>
    void setTheta(std::pair<II, II> thetas_) {
        std::copy(thetas_.first, thetas_.second, std::back_inserter(getModifyThetas()));
    }

    /**
     * @brief const getter
     */
    const std::vector<fvector_t>& getThetas() const {
        return thetas;
    }

    /**
     * @brief non-const getter.
     */
    std::vector<fvector_t>& getModifyThetas() {
        return thetas;
    }

};
} //ns
