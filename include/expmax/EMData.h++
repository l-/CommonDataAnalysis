/**
 * @file EMData.h++
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
 * @class EMData
 *
 * template param datapoint_t
 *
 * @brief Any concrete EM subclass "has_a" EMData of the appropriate type (i.e. mono- or multivariate)
 */
template<class T>
class EMData {

public:

    /**
     * double or fvector_t
     */
    typedef T datapoint_t;

protected:

    /**
     * @brief Data dimensionality
     */
    const int D;

    /**
     * Reimplemented for all univariate variants
     */
    std::vector<datapoint_t> data;

public:

    /**
     * @brief Constructor.
     *
     * param[in] D_ data dimensionality
     *
     */
    EMData(const unsigned D_ = 1);

    /**
     * @brief Const getter
     */
    const std::vector<datapoint_t>& getData() const {
        return data;
    }

    /**
     * @brief Get one data point
     *
     * @param[in] n which one (0<=n<getN())
     *
     * @return the data point
     */
    const datapoint_t& getData(const unsigned n) const {
        return data[n];
    }

    /**
     * @brief Const getter
     */
    std::vector<datapoint_t>& getData() {
        return data;
    }

    /**
     * @brief Get number of data points
     */
    size_t getN() const {
        return data.size();
    }

    /**
     * @brief Get D (or 1 for univariate instances)
     *
     * @return D
     */
    size_t getDataDimensionality() const { return D; }

    /**
     * @brief Call this from setData
     */
    template<class II>
    void setDataProper(std::pair<II, II> data_) {
        std::copy(data_.first, data_.second, std::back_inserter(data));
        std::cout << typeid(this).name() << ": accepted " << data.size() << " data points.\n";
    }

};

} // namespace
