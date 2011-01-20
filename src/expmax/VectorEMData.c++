/**
 * @file VectorEMData.c++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 20, 2011
 *
 */

#include "expmax/VectorEMData.h++"

using namespace CDA;

size_t VectorEMData::getDataDimensionality() const { return D; }

size_t VectorEMData::getD() const { return D; }

VectorEMData::VectorEMData(const VectorEMData& other)
  : EMData(other), D(other.D) {}

VectorEMData::VectorEMData(const unsigned D_) : D(D_) {};
