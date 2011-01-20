/**
 * @file EMGenericMixtureModelCore.c++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 */

#include "expmax/EMGenericMixtureModelCore.h++"

using namespace CDA;

template <class theta_T>
void EMGenericMixtureModelCore<theta_T>::update_thetas() {

#ifdef VERBOSE
    std::cout << className() << ": M step.\n";
#endif

    // The common part: updating of overall class probabilities.
    // It is assumed to be residing at [0] of the parameter vector.
    unsigned int K = getK();

    double sumclassif[K];
    for (unsigned k=0; k<K; ++k) {
        sumclassif[k] = 0.0; }

    for (unsigned n=0; n<getN(); ++n) {
        for (unsigned k=0; k<getK(); ++k) {
            sumclassif[k] += getClassif(n)(k);
        }
    }

    for (unsigned k=0; k<K; ++k) {
        getThetaObj().getModifyThetas(k,0) = 1/((double)(getN())) * sumclassif[k];

#ifdef VERBOSE
        std::cout << "Alpha " << k << " now " << getThetaObj().getClassProb() << std::endl;
#endif
    }

    improveClusterModelParameters();
}

template <class theta_T>
const std::string EMGenericMixtureModelCore<theta_T>::className() const {
    return std::string("EMGenericMixtureModelCore"); // shouldn't occur, class is still abstract
}

template <class theta_T>
const std::string EMGenericMixtureModelCore<theta_T>::getCSVHeader() const {
    std::stringstream out;

    out << "iteration" << ";";
    out << "k" << ";";
    out << "p" << ";";
    for (unsigned p=0; p<P; ++p) {
        out << paramName(p);
        if (p < P - 1) {
            out << ";";
        } else {
            out << "\n";
        }
    }

    return out.str();
}
