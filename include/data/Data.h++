/**
 * @file Data.h++
 * @version 0.02
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Dec 30, 2010
 *
 * @section LICENSE
 * LGPL v3.0
 *
 * @section COMMENT
 * Sense the POWER!
 *
 * Beginnings of a GNU R-like dataframe infrastructure.
 * Very unfinished
 */

#pragma once

#include <boost/foreach.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>
#include <iostream>
#include <algorithm>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/covariance.hpp>

#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/function.hpp>
// #include <boost/lambda/lambda.hpp>
// #include <boost/lambda/bind.hpp>
#include <boost/mem_fn.hpp>

#include <boost/bind.hpp>
#include <boost/bind/arg.hpp>
#include <boost/bind/placeholders.hpp>
#include <boost/bind/mem_fn.hpp>

#include "common_definitions.h++"

namespace CDA {

using namespace boost::accumulators;
namespace ublas = boost::numeric::ublas;

class EnhancedDatasetView {

};

/**
 * @class VectorComponentExtractor
 * @brief This is much clearer than the dopey boost::function version anyway
 */
struct VectorComponentExtractor {
  /**
   * @brief So as to conform with interface required by transform_iterator
   */
  typedef fvector_t::value_type result_type;
  const int i;
  VectorComponentExtractor(const int i_) : i(i_) {}
  fvector_t::value_type operator()(const fvector_t& v) { return v(i); }
};


/**
 * @class EnhancedDataset
 *
 * @brief aware of basic statistics for normalization etc.
 */
class MYEXPORT EnhancedDataset {

    // @todo Assign data types (for example, is it a real number, an angle, etc.)

public:
    typedef std::vector<fvector_t> dataframe_t; // @todo nummern auch bei filterung ... !
    typedef fvector_t row_t;
    dataframe_t m_data;
    typedef accumulator_set<double, stats<tag::mean, tag::variance, tag::max, tag::min> > boost_statistics_t;
    unsigned D;

    /**
     * @brief One for each dimension
     */
    std::vector<boost_statistics_t> accs;

    EnhancedDataset() : D(0) {}

    EnhancedDataset(unsigned D_)
    : D(D_) {
        for (unsigned d = 0; d < D; ++d) {
            boost_statistics_t x;
            accs . push_back( x );
        }
    }

    void addRow(const fvector_t& datarow) {
        assert(datarow.size() == D);

        m_data.push_back(datarow);

        for (unsigned i=0; i<D; ++i) {
            accs[i](datarow[i]);
        }
    }

    const dataframe_t& getData() const {
        return m_data;
    }

    size_t N() const { return m_data.size(); }

    template<class boost_accumtag>
    fvector_t getStatistics() const {
        fvector_t result(D);
        for (unsigned d = 0; d < D; ++d) {
            result(d) = extract_result< boost_accumtag >( accs[d] );
        }
        return result;
    }

    /**
     * @brief Handy for quick normalization, all dimensions scaled by same factor
     */
    double getMinOfAny() const {
        fvector_t presult = getStatistics<tag::min>();
        return *std::min_element(presult.begin(), presult.end());
    }

    /**
     * @brief Handy for quick normalization, all dimensions scaled by same factor
     */
    double getMaxOfAny() const {
        fvector_t presult = getStatistics<tag::max>();
        return *std::min_element(presult.begin(), presult.end());
    }

    fvector_t getNormalized(unsigned idx) const {
        return normalize(m_data[idx]);
    }

    const fvector_t& normalize(const fvector_t& in) const {
        fvector_t result(D);
        for (unsigned d = 0; d < D; ++d) {
            double max = extract_result< tag::max >( accs[d] );
            double min = extract_result< tag::min >( accs[d] );
            result(d) = (in(d) - min)/(max-min);
        }
        return result;
    }

    template<class indexcoll_t>
    fvector_t normalize_selective(const fvector_t& in, indexcoll_t& slc) const {

        // naja
        fvector_t result = in;

        BOOST_FOREACH (unsigned int d, slc) {
            double max = extract_result< tag::max >( accs[d] );
            double min = extract_result< tag::min >( accs[d] );
            result(d) = (in(d) - min)/(max-min);
        }

        return result;
    }

    /**
     * @brief Normalize only the selected fields
     */
    template<class indexcoll_t>
    boost::shared_ptr<EnhancedDataset> normalized(indexcoll_t& slc) const {

        boost::shared_ptr<EnhancedDataset> result( new EnhancedDataset(D) );

        for (dataframe_t::const_iterator dit = m_data.begin();
                dit != m_data.end();
                ++dit) {

            result -> addRow( normalize_selective(*dit, slc) );
        }

        assert(result->D == D);

        return result;
    }

    // Must be ref. because copy constructor neither implemented nor good
    template<class indexcoll_t>
    boost::shared_ptr<EnhancedDataset> projected(indexcoll_t& slc) const {
        // yes ... improve someday. do several in one go

        boost::shared_ptr<EnhancedDataset> result( new EnhancedDataset(slc.size()) );

        for (dataframe_t::const_iterator dit = m_data.begin();
                dit != m_data.end();
                ++dit) {

            fvector_t row(slc.size());

            int i=0;
            for (typename indexcoll_t::const_iterator slit = slc.begin();
                    slit != slc.end();
                    ++slit) {
                row(i++) = (*dit)(*slit);
            }

            result -> addRow(row);
        }

        return result;
    }

    typedef boost::function<fvector_t(const fvector_t&)> vector_transform_t;

    // vv mit ublas::slice
    //
    typedef boost::transform_iterator<vector_transform_t,
            dataframe_t::const_iterator
            > vtrit_t;
    //
    std::pair<vtrit_t, vtrit_t> getNormalizedIterPair() const {
        //
        //        // http://stackoverflow.com/questions/3413044/declaring-and-defining-a-function-object-inside-a-class-member-function

        boost::function<const fvector_t&(const fvector_t&)> fn
                (boost::bind(boost::mem_fn(&EnhancedDataset::normalize), this, _1 ));

        return std::make_pair(
                boost::make_transform_iterator(m_data.begin(), fn),
                boost::make_transform_iterator(m_data.end(), fn)
        );
    }

    /**
     * Transform data vectors by user-defined function!
     */
    std::pair<vtrit_t, vtrit_t> getTransformedIterPair(vector_transform_t& f) const {
        return std::make_pair(
                boost::make_transform_iterator(m_data.begin(), f),
                boost::make_transform_iterator(m_data.end(), f)
        );
    }

    typedef boost::transform_iterator<boost::function<fvector_t::value_type(const fvector_t&)>,
            dataframe_t::const_iterator
            > vtrit_1d_t;
    /**
     * @brief Extract one dimension only, in case the user needs a simple list of scalars
     *
     * @param[in] i select this dimension.
     */
    std::pair<vtrit_1d_t, vtrit_1d_t> getScalarIterPair(const int i) const {
        //
        //        // http://stackoverflow.com/questions/3413044/declaring-and-defining-a-function-object-inside-a-class-member-function

//        typedef const fvector_t::value_type& crap;
//        // what kind of unsightly stuff is this???
//        // What's more, it only JUST works with G++ and not at all with MSVC ...
//        typedef crap (fvector_t::*ttt)(unsigned int) const;
//
//        boost::function<crap(const fvector_t&)> fn
//                (boost::bind(
//                        boost::mem_fn(
//                             // This way it is not a "non-member function type"
//                             static_cast<ttt>(
//                             &fvector_t::operator())), _1, i ));
//                        // boost::mem_fn((const fvector_t::value_type&(const fvector_t&, unsigned int))(&fvector_t::operator())), _1, i ));

        VectorComponentExtractor fn(i);

        return std::make_pair(
                boost::make_transform_iterator(m_data.begin(), fn),
                boost::make_transform_iterator(m_data.end(), fn)
        );
    }


    /**
     * @brief Get a simple iterator pair
     *
     * @return
     */
    std::pair<dataframe_t::const_iterator, dataframe_t::const_iterator> getIterPair() const {
        return std::make_pair(m_data.begin(), m_data.end());
    }

    // EnhancedDatasetView projection()

    template<class F>
    std::pair<boost::filter_iterator<F, dataframe_t::const_iterator>,
              boost::filter_iterator<F, dataframe_t::const_iterator>
             > getFilteredIterPair(F& function) const {
        return std::make_pair(boost::make_filter_iterator(m_data.begin(), function),
                              boost::make_filter_iterator(m_data.end(), function));
    }

//    template <class Coll>
//    std::pair<boost::filter_iterator<boost::function<bool(dataframe_t::value_type)>, dataframe_t::const_iterator>,
//              boost::filter_iterator<boost::function<bool(dataframe_t::value_type)>, dataframe_t::const_iterator> >
//        getFilteredIterPair_above(Coll ) const {
//    // Coll is a  list of thresholds
//        boost::function<bool(const fvector_t&)> fn
//                        (boost::bind(std::less<double>,  boost::mem_fn(&EnhancedDataset::normalize), this, _1 ));
//
//    }
    //

    //    std::pair<boost::filter_iterator<boost::less>, dataframe_t::iterator>
    //      getFilteredIterPair_UpwardsFromNStddev() const {
    //
    //    }

private:
    EnhancedDataset(const EnhancedDataset&) {
        const bool no_copy_constructor_yet = false;
        assert(no_copy_constructor_yet);
    }
};

} // namespace
