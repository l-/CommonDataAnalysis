/**
 * @file VectorEMData.h++
 * @version 0.1
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 20, 2011
 *
 * when data comes in vectors of floating-point values
 */

#pragma once

#include "EMData.h++"

namespace CDA {

/**
 * @class VectorEMData
 *
 * @brief Data which comes in vectors.
 */
class VectorEMData : public EMData<fvector_t> {

protected:

    /**
     * @brief Yes?
     */
    using EMData<fvector_t>::data;

    /**
     * @brief Data dimensionality
     */
    const int D;

public:

    /**
     * @brief Obvious
     *
     * @param[in] D_ dimensionality
     */
    VectorEMData(const unsigned D_);

    /**
     * @brief Do the obvious
     */
    VectorEMData(const VectorEMData& other);

    /**
     * @brief Get D (or 1 for univariate instances)
     *
     * @return D
     */
    size_t getDataDimensionality() const;

    /**
    * @return D
    */
    size_t getD() const;

};

} // namespace
