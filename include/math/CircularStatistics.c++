/**
 * @file CircularStatistics.c++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 5, 2011
 *
 */

#include "CircularStatistics.h++"

using namespace CDA;

double wrap2pi(const double phi)
{
  return fmod(phi, 2*M_PI);
}
