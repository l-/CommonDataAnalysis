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

void EMGenericMixtureModelCore::update_thetas() {

#ifdef VERBOSE
    std::cout << className() << ": M step.\n";
#endif

    // The common part: updating of overall class probabilities
    double sumclassif[K];
    for (unsigned k=0; k<K; ++k) {
        sumclassif[k] = 0.0; }

    for (unsigned n=0; n<getN(); ++n) {
        for (unsigned k=0; k<getK(); ++k) {
            sumclassif[k] += classif[n](k);
        }
    }

    for (unsigned k=0; k<K; ++k) {
        m_theta.getModifyThetas()[k](0) = 1/((double)(getN())) * sumclassif[k];

#ifdef VERBOSE
        std::cout << "Alpha " << k << " now " << m_theta.getThetas()[k](0) << std::endl;
#endif
    }

    improveClusterModelParameters();
}

const std::string EMGenericMixtureModelCore::className() const {
    return std::string("EMGenericMixtureModelCore"); // shouldn't occur, class is still abstract
}


const std::string EMGenericMixtureModelCore::getCSVHeader() const {
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
