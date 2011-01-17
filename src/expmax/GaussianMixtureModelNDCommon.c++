/**
 * @file GaussianMixtureModelNDCommon.c++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 */

#include "expmax/GaussianMixtureModelNDCommon.h++"

using namespace CDA;

const std::string GaussianMixtureModelNDCommon::paramName(const unsigned p) const {
    std::stringstream out;
    if (p <= m_data . getDataDimensionality()) {
        out << "m_" << p;
    } else {
        int a = p - m_data . getDataDimensionality();
        out << "s_" << i(a) << "_" << j(a) << std::endl;
    }
    return out.str();
}
