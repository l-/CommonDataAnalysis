/**
 * @file ProbabilisticClustering.h++
 * @author  Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 * @version 0.1
 *
 * @section LICENSE
 *
 * This file is released under LGPL v3.0.
 *
 * @section DESCRIPTION
 *
 * Gaussian Mixture Model for fitting of a known number of Gaussians;
 * and other probabilistic models of use in the processing of images.
 *
 * Common patterns are extracted afap to build a hierarchy of classes.
 *
 */

#pragma once

#include <typeinfo> // Just to print names of classes, I promise

#include <boost/foreach.hpp>
#include <boost/optional.hpp>

// Type assertions.
#include <boost/mpl/equal.hpp>

#include <boost/function.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
using namespace boost::lambda;

// Loop-avoiding trickery
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/algorithm/minmax_element.hpp>

// Sine qua non of modern C++
#include <boost/shared_ptr.hpp>

// #include <boost/array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

// Numeric limits
#include <limits>

#include <iostream> // argh, verflechtung

#include "math/vonMises.h++"
#include "math/numeric_optimization.h++"
#include "utils/comparisons.h++"

/**
 * @brief Common declaration: Vector type.
 */
typedef boost::numeric::ublas::vector<double> fvector_t;

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

/**
 * @class EMThetas
 *
 * @brief Model parameters for EM
 *
 * @section DESCRIPTION
 * To avoid duplication. Can't just put it into EM and inherit
 * because of duplicate class problems
 *
 * @section ARCHITECTURE CAVEAT
 * I dislike having this class around, this should all be in EM
 * in any proper programming language, it would be possible.
 */
class EMThetas {

protected:
    /**
     * @brief Model parameters
     * possibly a single vector, or K PDF parameter vectors, one per class
     * => these belong in the concrete classes? I hate duplication!
     */
    std::vector<fvector_t> thetas;

public:
    /**
     * Initial guess for model parameters -- strictly needed
     *
     * Input data format differs for each individual distribution,
     * be careful and read the documentation.
     *
     * @todo random numbers when none given
     */
    template<class II>
    void setTheta(std::pair<II, II> thetas_) {
        std::copy(thetas_.first, thetas_.second, std::back_inserter(getModifyThetas()));
    }

    /**
     * @brief const getter
     */
    const std::vector<fvector_t>& getThetas() const {
        return thetas;
    }

    /**
     * @brief non-const getter.
     */
    std::vector<fvector_t>& getModifyThetas() {
        return thetas;
    }

};


/**
 * @class EM
 *
 * @brief A generic EM model fitting class.
 * Contains the general outline of EM algorithm,
 * details to be filled in by implementations (not that easy in C++ ...)
 */
template<class datapoint_t>
class EM {

protected:

    /**
     * @brief In C++, _always_ works better than inheritance. I should have remembered ...
     */
    EMThetas m_theta;

    /**
     * @brief In C++, _always_ works better than inheritance
     */
    EMData<datapoint_t> m_data;

    /**
     * @brief One iteration.
     */
    double EMstep() {
        update_hidden();
        update_thetas();
        return logLikelihood();
    }

    /**
     * @brief <b>E-step</b>: update hidden attributes, i.e. class membership probabilities
     */
    virtual void update_hidden() = 0;

    /**
     * @brief <b>M-step</b>: update parameters of model PDF
     */
    virtual void update_thetas() = 0;

    /**
     * @brief Important. This is being optimized
     */
    virtual double logLikelihood() const = 0;

    /**
     * @brief Implementer should overwrite this method, since typeid(this).name() returns crap
     */
    virtual const std::string className() const {
        return typeid(this).name();
    }

public:

    /**
     * @brief Constructor, to be called explicitly in multivariate case
     */
    EM(const unsigned D_ = 1)
    : m_data(D_) { }

    /**
     * @brief Getter for estimate
     */
    virtual const fvector_t& getHiddenParamEstimate(const unsigned n) const = 0;

    /**
     * @brief Initial guess for model parameters -- strictly needed
     *
     * Input data format differs for each individual distribution,
     * be careful and read the documentation.
     *
     * @param[in] thetas
     */
    template<class II>
    void setTheta(std::pair<II, II> thetas) {
        m_theta . setTheta(thetas);
    }


    /**
     * @brief Run Expectation Maximization
     *
     * @param[in] MAXITER up to a certain number of iterations
     * @param[in] thresh  stop when loglikelihood difference below threshold
     * @param[in] output_csv
     *
     */
    void EMrun(const unsigned MAXITER, const double thresh, const boost::optional<std::ostream*> output_csv = boost::none);

    /**
     * @brief Output CSV
     *
     * @param[in] iteration (optional) output iteration no. in first column
     *
     */
    virtual const std::string dumpParameters(const boost::optional<unsigned int> iteration) const = 0;

    /**
     * @brief Möchte es nicht zu kompliziert machen (14/12/2010),
     * daher bekommen die ganzen Klassen erst einmal CSV-Ausgabe hardwired.
     */
    virtual const std::string getCSVHeader() const = 0;

    /**
     * @brief Initialize known samples, i.e. your data points
     * Also initialize the classif vector
     */
    template<class II>
    void setData(std::pair<II, II> data_) {

        // II must iterate over datapoint_t elements
        BOOST_MPL_ASSERT_MSG((boost::mpl::equal<datapoint_t, typename II::value_type>::type::value), UnsupportedDatavectorType, (typename II::value_type));

        m_data . setDataProper(data_);
    }
};

/**
 * @class FitMulticlassByEM
 *
 * An abstract model fitting class
 *
 * template param: the datapoint_t
 *
 */
template<class T>
class FitMulticlassByEM : public EM<T> {

public:
    /**
     * @brief Propagate the datapoint_t.
     * Why can't I give it the same name? Annoying.
     */
    typedef T datapoint_t;

    using EM<datapoint_t>::className;

protected:

    /**
     * @brief Number of individual Distributions
     */
    const unsigned K;

    /**
     * @brief N stück: class appartenance probabilities
     */
    std::vector<fvector_t> classif;

    using EM<datapoint_t>::m_theta;
    using EM<datapoint_t>::m_data;

    /**
     * @brief <b>E-step</b>: update hidden attributes, i.e. classification
     *
     * @section NOTE
     * It is common to all multiclass (mixture model) problem formulations.
     */
    void update_hidden();

public:

    /**
     * Constructor
     *
     * To be called by subclasses
     *
     * @param[in] D_ Data Dimensionality
     * @param[in] K_ Number of Classes (~)
     *
     */
    FitMulticlassByEM(const unsigned K_, const unsigned D_ = 1)
    : EM<datapoint_t>(D_)
      , K(K_)
      { }

    /**
     * @brief Output CSV
     *
     * @param[in] k class
     * @param[in] iteration (optional) output iteration no. in first column
     *
     * @return one CSV line
     */
    const std::string dumpParameters(const unsigned k, const boost::optional<unsigned int> iteration) const;

    /**
     * @brief Output CSV
     *
     * @param[in] iteration (optional) output iteration no. in first column
     *
     * @return multiple CSV lines
     */
    const std::string dumpParameters(const boost::optional<unsigned int> iteration) const;

    /**
     * @brief Conveniently encapsulate access ofm_theta
     *
     * all parameters except the outer probability
     */
    double getParam(const unsigned k, const unsigned p) const;

    /**
     * @brief Header to go with the CSV output. Must be reimplemented
     */
    virtual const std::string getCSVHeader() const = 0;

    /**
     * @brief Implemented here.
     */
    double logLikelihood() const;

    /**
     * @brief Getter for estimate
     *
     * @param[in] n Point number n
     */
    const fvector_t& getHiddenParamEstimate(const unsigned n) const {
        return classif[n];
    }

    /**
     * @brief Better have a getter for this too
     */
    unsigned int getK() const {
        return K;
    }

    /**
     * @brief Get class with best estimated membership probability
     *
     * @param[in] n Point number n
     */
    unsigned int getBestClass(const unsigned n) const {
        const fvector_t& cs = getHiddenParamEstimate(n);
        return *( std::max_element<boost::counting_iterator<int>, compareByVectorElementValue<fvector_t, size_t> >
        (boost::counting_iterator<int>(0),
                boost::counting_iterator<int>(getK()),
                compareByVectorElementValue<fvector_t>(cs)) );
    }

    /**
     * @brief Read classification results (1)
     *
     * @return Best class for each point, in the same order
     * as the points were entered (so you can just access 'em by index,
     * as everywhere)
     */
    std::vector<unsigned int> getClassifList() const {
        std::vector<unsigned int> result;
        std::for_each(boost::counting_iterator<int>(0),
                boost::counting_iterator<int>(m_data.getN()),
                boost::lambda::bind(boost::mem_fn(&std::vector<unsigned int>::push_back),
                        &result,
                        boost::lambda::bind(
                                boost::mem_fn(&FitMulticlassByEM<datapoint_t>::getBestClass),
                                this,
                                boost::lambda::_1)));
        return result;
    }

    /**
     * @brief Get current estimated class weight
     *
     * @param[in] k class no.
     */
    double getPk(const unsigned k) const;

    /**
     * @brief Initialize known samples, i.e. your data points
     * Also initialize the classif vector
     */
    template<class II>
    void setData(std::pair<II, II> data_) {

        // II must iterate over datapoint_t elements
        BOOST_MPL_ASSERT_MSG((boost::mpl::equal<datapoint_t, typename II::value_type>::type::value), UnsupportedDatavectorType, (typename II::value_type));

        m_data . setDataProper(data_);
        initClassif();
    }

    /**
     * @brief Initialize class membership beliefs
     *
     * @section NOTE
     * setData of derived classes must call this
     * => done, defined setData here. can one be "using" a template fn?
     */
    void initClassif() {
        classif . clear();

        for (unsigned i=0; i<m_data.getN(); ++i) {
            classif . push_back ( fvector_t(K, 1/(double)K) );
        }
    }

    /**
     * @brief Evaluate the cluster model PDF with class parameters \f$k\f$ at \f$\vec{x}\f$
     *
     * @param[in] k class
     * @param[in] x feature vector
     *
     * @return probability \f$p(x,k\vert\theta)\f$
     */
    virtual double evalPDF(const unsigned k, const datapoint_t& x) const = 0;

};

/**
 * @brief It's been demoted to a typedef.
 */
typedef FitMulticlassByEM<double> FitUnivariateMulticlassByEM;

/**
 * @brief It's been demoted to a typedef.
 */
class FitMultivariateMulticlassByEM : public FitMulticlassByEM<fvector_t> {

public:

    /**
     * @brief Define the datapoint_t
     */
    typedef fvector_t datapoint_t;

public:

    /**
     * @brief Constructor
     *
     * @param[in] K_ this way round, because of default parameter in superclass
     * @param[in] D_
     */
    FitMultivariateMulticlassByEM(const unsigned K_, const unsigned D_)
    : FitMulticlassByEM<fvector_t>(K_, D_) {}
    // , EM<fvector_t>(D_) {}

    /**
     * @brief Convenience getter for D, same format as getN() and getP()
     */
    unsigned getD() const;
};

/**
 * @class GaussianMixtureModel1D
 *
 * A GaussianMixtureModel class, with EM fitting routine built in
 *
 */
class GaussianMixtureModel1D : public FitUnivariateMulticlassByEM {

public:
    typedef double datapoint_t;

private:
    using EM<datapoint_t>::m_theta;
    using EM<datapoint_t>::m_data;

    using FitUnivariateMulticlassByEM::K;

    inline double squaredDistanceToMean(const unsigned k, const datapoint_t x) const {
        return fabs(x - getMean(k)); // OBS! fabs!
    }

protected:

    /**
     * @brief <b>E-step</b>: update hidden attributes, i.e. class membership probabilities
     */
    void update_thetas();

public:

    /**
     * @brief Evaluate the cluster model PDF with class parameters \f$k\f$ at \f$\vec{x}\f$
     *
     * @param[in] k class
     * @param[in] x feature vector
     *
     * @return probability \f$p(x,k\vert\theta)\f$
     */
    inline double evalPDF(const unsigned k, const datapoint_t& x) const;

    /**
     * @brief Call superclass constructor to fill data fields
     */
    GaussianMixtureModel1D(const unsigned K_)
    : FitUnivariateMulticlassByEM(K_)
    {  }

    /**
     * @brief Estimated parameter getter.
     * Get current estimated mean vector of class k
     *
     * @param[in] k class no.
     */
    double getMean(const unsigned k) const;

    /**
     * @brief Estimated parameter getter
     * Get current estimated sigma of class k
     *
     * @param[in] k class no.
     */
    double getSigma(const unsigned k) const;
    //
    //    /**
    //     * Evaluate the PDF with parameters of class k
    //     *
    //     * @param[in] k class no.
    //     * @param[in] x
    //     */
    //    inline double evalPDF(const unsigned k, const double x) const;

    /**
     * @brief Initialize known samples, i.e. your data points
     *
     * @param[in] data_ A pair of iterators.
     * @todo we should check II for compliance.
     */
    template<class II>
    void setData(std::pair<II, II> data_) {
        FitUnivariateMulticlassByEM::setData(data_);
    }

    /**
     * @brief Find discriminant between two classes
     *
     * @param[in] k1 a class
     * @param[in] k2 a class
     *
     * At least this. Is there an easier way?
     */
    boost::function<double(const double)>
    getDiscriminantFunction(const int k1, const int k2) const {
        // @todo double-check

        assert(k1 != k2);

        return boost::function<double(const double)>(
                boost::lambda::bind(fabs, (boost::lambda::_1 - getMean(k1))) * pow(getSigma(k2),2)
                - boost::lambda::bind(fabs, (boost::lambda::_1 - getMean(k2)))  * pow(getSigma(k1),2)
                - 2 * log ( getSigma(k2) / getSigma(k1) ) * pow(getSigma(k2),2) * pow(getSigma(k1),2));
        //          < 0  ? k1 : k2 );

        // eh? ist auch falsch ...
    }

    /**
     * @brief Find explicit class boundary
     *
     * @param[in] k1 Class 1
     * @param[in] k2 Class 2
     *
     * @return x>value: in class 2
     */
    double findDecisionLimit(const int k1, const int k2) const {

        double min = getMean(k1);
        double max = getMean(k2);

        if (max > min)
            return
            findSingleUnivariateRootIntervalSearch(getDiscriminantFunction(k1, k2), min, max,1e-7);
        else
            return
            findSingleUnivariateRootIntervalSearch(getDiscriminantFunction(k2, k1), max, min, 1e-7);

        // @muß getestet werden
    }

    /**
     * @brief Class name for logging ...
     */
    const std::string className() const;

    /**
     * @brief Each class gets its own line, that's the easiest way
     */
    const std::string getCSVHeader() const;

};

/**
 * @class EMGenericMixtureModelCore
 *
 * @brief Ideally, each model knows its data layout, parameter estimators and stuff.
 * There could be a general blueprint for dealing with models which cannot be handled
 * analytically, prompting a gradient search for the MLE step
 *
 * @section NOTE
 * for multivariate.
 */
class EMGenericMixtureModelCore : public FitMultivariateMulticlassByEM {

public:
    /**
     * @brief Propagate the datapoint_t
     */
    typedef FitMultivariateMulticlassByEM::datapoint_t datapoint_t;

protected:

    /**
     * @brief Parameter space dimensionality.
     *
     * @section NOTA BENE
     * We don't count the class probability here,
     * as it is handled differently.
     */
    const unsigned P;

    /**
     * @brief <b>M-step</b>: update parameters of model PDF
     *
     * It calls improveClusterModelParameters, where the implementer
     * specifies how the other parameters (beside a priori class probabilities)
     * are updated.
     */
    void update_thetas();

    /**
     * @brief The quasi-MLE step. Since this may be done in different-different ways,
     * it can be mixed in from a separate class.
     *
     * It is called by update_thetas
     */
    virtual void improveClusterModelParameters() = 0;

    /**
     * @brief Get param space dimensionality.
     *
     * @return If it is still zero, there's a programming error!
     */
    unsigned int getP() const { assert(P!=0); return P; }

public:

    /**
     * @brief Constructor to be called by implementer, please
     *
     * @param[in] K_
     * @param[in] D_ this way round because of default param further up
     * @param[in] P_ specific to the model.
     */
    EMGenericMixtureModelCore(const unsigned K_, const unsigned D_, const unsigned P_)
    : FitMultivariateMulticlassByEM(K_, D_), P(P_) {}

    /**
     * @brief The individual way of evaluating the PDF
     *
     * @param k class no.
     * @param x
     *
     * @return probability \f$p(x,k\vert\theta)\f$
     */
    virtual double evalPDF(const unsigned k, const datapoint_t& x) const = 0;

    /**
     * @brief Class name for logging ...
     */
    virtual const std::string className() const;

    /**
     * @brief Name of a parameter p=\f$\theta_p\f$
     * (For output)
     *
     * @param[in] p 0<=p<P
     */
    virtual const std::string paramName(const unsigned p) const = 0;

    /**
     * @brief Each class gets its own line, that's the easiest way
     * Assume iteration no. is to be shown always
     *
     * @return String the parameter names together.
     */
    const std::string getCSVHeader() const;

};

/**
 * @class EMMbyGradientDescent
 *
 * @brief EM whose M-step is done by gradient descent
 */
class EMMbyGradientDescent : public EMGenericMixtureModelCore {

protected:

    /**
     * @brief Call superclass constructor only.
     */
    EMMbyGradientDescent(const unsigned D_, const unsigned K_, const unsigned P_)
    : EMGenericMixtureModelCore(D_, K_, P_) {}

    /**
     * @brief The quasi-MLE step (the not-so-easy part)
     *
     * Implemented here.
     */
    void improveClusterModelParameters();

    /**
     * @brief if there are any values which are needed later
     * @todo refactor!
     */
    virtual void preparations_for_evalPDFderivP() = 0;

    virtual double evalPDFderivP(const unsigned k, const datapoint_t& x, const unsigned p) const = 0;

    /**
     * @brief Now if we are able to evaluate the PDF, even better: calculate a gradient,
     * we can actually do something.
     *
     * @section NOTE
     * This is used in the M-step to improve on all parameters.
     *
     * @param[in] p 0<=p<P
     *
     * @return The component of the gradient \f$\frac{\partial L}{\partial p}\f$ <b>(nicht wahr?)</b>
     *
     */
    double improvementOnTheta(const unsigned p) const;
};

/**
 * @class GaussianMixtureModel
 *
 * @brief A heteroscedastic N-dimensional Gaussian Mixture Model
 */
class GaussianMixtureModel : public EMMbyGradientDescent {

private:

    typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> sym_mtx_t;
    std::vector<sym_mtx_t> m_cached_invsigmas;

    /**
     * @brief @todo <b>This is still all pretty naïve, unsophisticated, slow ...</b>
     * But I've no nerve for that kind of optimization right now
     */
    double getSigmaCofactor(const unsigned k, const unsigned i, const unsigned j) const;

    /**
     * @brief @todo <b>This is still all pretty naïve, unsophisticated, slow ...</b>
     * But I've no nerve for that kind of optimization right now
     */
    double getInvSigmaDerivTerm(const unsigned k, const unsigned i, const unsigned j, const datapoint_t&) const;

    /**
     * @brief at least some optimization, otherwise it is ridiculous
     */
    void preparations_for_evalPDFderivP();

    /**
     * @brief Get param index - D of covariance parameter sigma(i,j) ;-)
     *
     * @section NOTA BENE
     * <b> \f$i < j\f$ always! </b>
     * So the upper triangular part of the matrix is stored.
     */
    inline unsigned a(const unsigned i, const unsigned j) const {
        return i*m_data . getDataDimensionality() - j - (i*(i+1))/2;
    }

    /**
     * @brief Get i index of covariance parameter ;-)
     */
    inline unsigned i(const unsigned a, const unsigned iter = 1) const {
        const unsigned D = m_data . getDataDimensionality();
        if (a>=D) { return i(a-D+iter, iter+1); }
        else { return a; }
    }

    /**
     * @brief Get j index of covariance parameter ;-)
     */
    inline unsigned j(const unsigned a, const unsigned iter = 1) const {
        const unsigned D = m_data . getDataDimensionality();
        if (a>=D) { return j(a-D+iter, iter+1); }
        else { return iter-1; }
    }

public:

    /**
     * @brief Constructor
     *
     * @param[in] D_
     * @param[in] K_
     *
     * @section Parameter space dimensionality
     * P = D + D(D+1)/2 (mean + covariances)
     *
     */
    GaussianMixtureModel(const unsigned D_, const unsigned K_)
    : EMMbyGradientDescent(D_, K_, D_ + D_*(D_+1)/2) {}

    /**
     * @brief Construct names of parameters, e.g. for output
     */
    const std::string paramName(const unsigned p) const;

    /**
     * @brief Evaluate model PDF of cluster k
     *
     * @param k class no.
     * @param x feature vector \f$\vec{x}\f$
     *
     * @return \f$p(x,k\vert\theta) = \frac{1}{(2\pi)^{\frac{D}{2}}\left|\Sigma_k\right|}e^{\frac{-1}{2}(\vec{x}-\vec{\mu_k})^T\Sigma_k^{-1}(\vec{x}-\vec{\mu_k})}\f$
     *
     */
    double evalPDF(const unsigned k, const datapoint_t& x) const;

    /**
     * @brief Get numerical value of derivateive
     *
     * @param k class no.
     * @param x feature vector \f$\vec{x}\f$
     * @param p number of parameter (first D are mean, next D(D+1)/2 are covariances)
     *
     */
    double evalPDFderivP(const unsigned k, const datapoint_t& x, const unsigned p) const;

    /**
     * @brief Extract the covariance matrix of Gaussian no. k
     *
     * @param[in] k class no.
     *
     * @return a ublas matrix
     */
    const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> getSigmaMatrix(const unsigned k) const;

    /**
     * @brief
     *
     * @return det(S)
     *
     * @todo urgent avoid re-calculations.
     * I'm not saying "memoize" but "be watchful".
     */
    double getSigmaDet(const unsigned k) const;

    /**
     * @brief Get inverse of Sigma of Gaussian no. k
     */
    const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> getInvSigma(const unsigned k) const;

    /**
     * @brief Get Mean vector of Gaussian no. k
     */
    const boost::numeric::ublas::vector_range<const fvector_t> getMean(const unsigned k) const;

};


/**
 * @class CircularMixtureModel1D
 *
 * A univariate VonMises-MixtureModel class, with EM fitting routine built in
 *
 */
class CircularMixtureModel1D : public FitUnivariateMulticlassByEM {

    typedef double datapoint_t;

    using EM<datapoint_t>::m_theta;
    using EM<datapoint_t>::m_data;

    using FitUnivariateMulticlassByEM::K;

protected:

    /**
     * @brief <b>E-step</b>: update hidden attributes, i.e. class membership probabilities
     */
    void update_thetas();

public:

    /**
     * @brief Evaluate the cluster model PDF with class parameters \f$k\f$ at \f$\vec{x}\f$
     *
     * @param[in] k class
     * @param[in] x feature vector
     *
     * @return probability \f$p(x,k\vert\theta)\f$
     */
    double evalPDF(const unsigned k, const datapoint_t& x) const;

    /**
     * @brief Call superclass constructor to fill data fields
     */
    CircularMixtureModel1D(const unsigned K_)
    : FitUnivariateMulticlassByEM(K_)
    {  }

    /**
     * @brief Estimated parameter getter.
     * Get current estimated mean vector of class k
     *
     * @param[in] k class no.
     */
    double getMu(const unsigned k) const;

    /**
     * @brief Estimated parameter getter
     * Get current estimated sigma of class k
     *
     * @param[in] k class no.
     */
    double getKappa(const unsigned k) const;

    /**
     * Evaluate the PDF with parameters of class k
     *
     * @param[in] k class no.
     * @param[in] x
     */
    inline double evalPDF(const unsigned k, const double x) const;

    /**
     * @brief Initialize known samples, i.e. your data points
     */
    template<class II>
    void setData(std::pair<II, II> data_) {
        FitUnivariateMulticlassByEM::setData(data_);
    }

    /**
     * @brief Class name for logging ...
     */
    const std::string className() const;

    /**
     * @brief Each class gets its own line, that's the easiest way
     */
    const std::string getCSVHeader() const;

};


/* UTILITIES: */

/**
 * @brief
 * Use univariate mixture model to determine a threshold
 * The classical application, from BV1 lecture actually!
 */
template<class Iter>
boost::shared_ptr<GaussianMixtureModel1D> thresholdFinder(const std::pair<Iter, Iter> input,
        const boost::optional<std::ostream*> output_csv = boost::none
) {

    std::vector<fvector_t> init_theta_bimodal;

    fvector_t hintergrund(3);
    fvector_t vordergrund(3);

    std::pair<Iter, Iter> minmax =
            boost::minmax_element(input.first, input.second);

    // Init prob.
    hintergrund(0) = 0.5; // if you know differently, please initialize it accordingly
    vordergrund(0) = 1 - hintergrund(0);

    // Init mean
    hintergrund(1) = * ( minmax.first );
    vordergrund(1) = * ( minmax.second );

    // Init stddev
    hintergrund(2) = 0.5 * (*minmax.second - *minmax.first);
    vordergrund(2) = hintergrund(2);

    init_theta_bimodal.push_back(hintergrund);
    init_theta_bimodal.push_back(vordergrund);

    boost::shared_ptr<GaussianMixtureModel1D> bimodal (new GaussianMixtureModel1D(2) );

    bimodal -> setData(input);
    bimodal -> setTheta(std::make_pair(init_theta_bimodal.begin(), init_theta_bimodal.end()));
    bimodal -> EMrun(50, 0.0000001, output_csv); // @todo time

    return bimodal;
}


/**
 * @brief Calculate the determinant of a small matrix.
 */
template<class matrix_t>
double det(matrix_t& m) {

    if (m.size1() == 1) {
        return m(0,0);
    }

    namespace ublas = boost::numeric::ublas;

    typedef ublas::permutation_matrix<std::size_t> permatrix;

    ublas::matrix<double> A = m;
    permatrix pm(A.size1());
    int lures = ublas::lu_factorize(A, pm);

    if( lures != 0 ) {
        std::cerr << "LU factorization of " << m << " failed." << std::endl;
        std::cerr << "you are stuck with: " << A << pm << std::endl;
    }

    double res = 1.0;
    for (unsigned i=0; i<A.size1(); ++i) {
        res *= A(i,i);
    }
    return res;
}


///**
// * Note: Circular and other manifold statistics: Maybe faking them via Gausses
// * (and specially defined distances) is en
// */
