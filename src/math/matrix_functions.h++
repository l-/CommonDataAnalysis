/**
 * @file matrix_functions.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 17, 2011
 *
 */

#pragma once

namespace CDA {

/**
 * @brief Calculate the determinant of a small matrix.
 *
 * used internally
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

#ifdef DETAIL_VERBOSE_2
    std::cerr << "Determinant = " << res << std::endl;
#endif

    return res;
}

} // namespace
