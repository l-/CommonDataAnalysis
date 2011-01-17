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

#include <typeinfo>

namespace CDA {

/**
 * @class EMData
 *
 * template param datapoint_t
 *
 * @brief Any concrete EM subclass "has_a" EMData of the appropriate type (i.e. mono- or multivariate)
 * It encapsulates access to the original data points.
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
     * Constructor, to be called explicitly in multivariate case only
     *
     * @param[in] D_ data dimensionality
     */
    EMData(const unsigned D_ = 1);

    /**
     * @brief Might be redefined by a subclass in case it wants its own version of EMData ...
     * so be careful, ONLY access it via accessor
     *
     * @return m_data
     */
    EMData<datapoint_t>& getDataObj();

    /**
     * @brief Same in const
     * @return reference to datapoints-holding object
     */
    const EMData<datapoint_t>& getDataObj() const;

    /**
     * @brief Const getter
     */
    const std::vector<datapoint_t>& getData() const;

    /**
     * @brief Get one data point
     *
     * @param[in] n which one (0<=n<getN())
     *
     * @return the data point
     */
    const datapoint_t& getData(const unsigned n) const;

    /**
     * @brief Const getter
     */
    std::vector<datapoint_t>& getData();

    /**
     * @brief Get number of data points
     *
     * formerly getN(); this is now a function in the class possessing this one
     */
    size_t getNumberOfDataPoints() const;

    /**
     * @brief Get D (or 1 for univariate instances)
     *
     * @return D
     */
    size_t getDataDimensionality() const;

    /**
     * @brief Call this from setData
     *
     * @param[in] iterators on data entries to copy
     */
    template<class II>
    void setDataProper(std::pair<II, II> data_) {
        std::copy(data_.first, data_.second, std::back_inserter(data));
#ifdef VERBOSE
        std::cout << typeid(this).name() << ": accepted " << data.size() << " data points.\n";
#endif
    }

};

} // namespace
