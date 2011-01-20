/**
* @file EMData.c++
*
* @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
*
* Created on: Jan 17, 2011
*
*/

#include "expmax/EMData.h++"

using namespace CDA;

namespace CDA{
template<> EMData<double>::EMData(const unsigned D_) {}; // all's well

template<class T>
EMData<T>::EMData(const EMData<T>& other)
  // : D(other.D)
{
#ifdef EXTRA_VERBOSE
    std::cerr << "EMData copy constructor called!\n";
#endif
    std::copy(other.data.begin(), other.data.end(), std::back_inserter(data));
}


} // namespace

template <class datapoint_t>
const datapoint_t& EMData<datapoint_t>::getData(const unsigned n) const {
    return data[n];
}

template <class datapoint_t>
const std::vector<datapoint_t>& EMData<datapoint_t>::getData() const {
    return data;
}

template <class datapoint_t>
std::vector<datapoint_t>& EMData<datapoint_t>::getData() {
    return data;
}

template <class datapoint_t>
size_t EMData<datapoint_t>::getNumberOfDataPoints() const {
    return data.size();
}

template <class datapoint_t>
void EMData<datapoint_t>::clear() {
    data.clear();
}

// Again, without these two declarations, the file would be useless ;-)
template class EMData<double>;
template class EMData<fvector_t>;
