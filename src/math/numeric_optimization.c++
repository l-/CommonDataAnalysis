/**
 * @file numeric_optimization.c++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Dec 28, 2010
 *
 */

#include "numeric_optimization.h++"

using namespace CDA;

double findSingleUnivariateRootIntervalSearch(const boost::function<double(double)>& f,
                                              const double xa,
                                              const double xb,
                                              const double epsilon) {

    // @todo evaluate less, this isn't haskell => OK

    double upperbound = xb, lowerbound = xa;
    double fxupper = f(upperbound), fxlower = f(lowerbound);

    assert(fxupper >= 0 && fxlower <= 0 && !(fxupper==fxlower));

    if (fxlower == 0) { return lowerbound; }
    if (fxupper == 0) { return upperbound; }

    double bisect, newy;

    while (fxupper - fxlower > epsilon) {
        bisect = lowerbound + 0.5 * (upperbound - lowerbound);
        newy = f(bisect);

        if (newy == 0) { return newy; } // if, by chance ...

        if (newy < 0) {
            lowerbound = bisect;
            fxlower = newy;
        } else {
            upperbound = bisect;
            fxupper = newy;
        }
    }

    return (upperbound + lowerbound)/2;
}
