/**
 * @file common_definitions.h++
 *
 * @author Erik Flick <erik.flick [AETT] informatik.uni-hamburg.de>
 *
 *  Created on: Jan 5, 2011
 *
 */

#pragma once

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>

namespace CDA {

/**
 * @brief Common declaration: Vector type.
 */
typedef boost::numeric::ublas::vector<double> fvector_t;
typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> sym_mtx_t;

}
	
#ifdef _WIN32
#define M_PI 2 * acos(0.0)
#define isnan(x) ((x)!=(x))
// AAAAAAARGH! M$ hell ... #define INFINITY (+1.0)/(+0.0)
// const double temp = 1.0;
// const double INFINITY = temp/(temp-1.0);
// => 1>Projekt : error PRJ0002 : Fehler "1073807364" wurde von "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin\x86_amd64\cl.exe" zurckgegeben.
// (we are being user friendly here) => using official limits now.
#endif

#ifdef _WIN32
#define MYEXPORT __declspec(dllexport)
#else
#define MYEXPORT  
#endif