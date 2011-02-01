/**
 * @file EMGenericMixtureModelCore.c++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 */

#include "expmax/EMGenericMixtureModelCore.h++"
#include "expmax/GaussianMixtureModelNDParams.h++"
#include "expmax/EMTheta.h++"

using namespace CDA;

template <class theta_T>
void EMGenericMixtureModelCore<theta_T>::update_thetas() {

#ifdef VERBOSE
    std::cout << className() << ": MixtureModel M step.\n";
#endif

    // The common part: updating of overall class probabilities.
    // It is assumed to be residing at [0] of the parameter vector.
    unsigned int K = getK();
    std::vector<double> sumclassif(K); // MSVC

    for (unsigned k=0; k<K; ++k) {
        sumclassif[k] = 0.0; }

    for (unsigned n=0; n<getN(); ++n) {
        for (unsigned k=0; k<getK(); ++k) {
            sumclassif[k] += getClassif(n)(k);
        }
    }

    for (unsigned k=0; k<K; ++k) {
        getThetaObj() -> getModifyThetas(k,0) = 1/((double)(getN())) * sumclassif[k];

#ifdef VERBOSE_2
        std::cout << "Alpha " << k << "/" << getK() << " now " << getThetaObj() -> getThetas(k,0) << std::endl;

        std::cout << "Sumclassif " << k << " = " << sumclassif[k] << std::endl;
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
    // out << "p" << ";";
    for (unsigned p=0; p<getP(); ++p) {
        out << paramName(p);
        if (p < getP()-1) {
            out << ";";
        } else {
            out << "\n";
        }
    }

    return out.str();
}

// The 2 specializations
namespace CDA {
template
class EMGenericMixtureModelCore<GaussianMixtureModelNDParams>;
template
class EMGenericMixtureModelCore<EMThetas>;
} // namespace
