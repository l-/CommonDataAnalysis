/**
 * @file comparisons.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Dec 29, 2010
 *
 */

#pragma once

namespace CDA {

/**
 * @class compareByVectorElementValue
 *
 * @brief Use this as a <b>functor</b> to compare numbers a and b
 * by the values v[a] and v[b]
 *
 * A "strict weak ordering", as per STL
 */
template<class vector_t, class idx_t = size_t, class LT = std::less<typename vector_t::value_type> >
struct compareByVectorElementValue {

    /**
     * @brief *reference* to the vector
     */
    const vector_t& m_v;

    /**
     * @brief Keep a reference of the vector, use it locally.
     */
    compareByVectorElementValue(const vector_t& v)
    : m_v(v) {

    }

    /**
     * @brief What do you expect?
     */
    bool operator()(const idx_t a, const idx_t b) {
        return LT()(m_v(a), m_v(b));
    }

};

} // namespace
