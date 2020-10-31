#ifndef MOLBIOLIB_DOUBLEVALUEFUNCTIONS_H
#define MOLBIOLIB_DOUBLEVALUEFUNCTIONS_H 1

//////////////////////////////////////////////////////////////////////////////
// MolBioLib: A C++11 framework for rapidly developing bioinformatics tasks //
// Copyright (C)  2012  Massachusetts General Hospital                      //
//                                                                          //
// This program is free software: you can redistribute it and/or modify     //
// it under the terms of the GNU General Public License as published by     //
// the Free Software Foundation, version 3 of the License.                  //
//                                                                          //
// This program is distributed AS IS.                                       //
//////////////////////////////////////////////////////////////////////////////


#include "src/include/PrimitiveTypes.hpp"


/** \file DoubleValueFunctions.hpp
 * Output functions on a double value.
 * \fn double factorial(double n)
 * Find the minimum element.  Usage:
 * \code
 * double theFactorial = factorial(5);
 * \endcode
 * We output to double since this uses the gamma function.
 */
double factorial(double n)
{
    return tgamma(n + 1.0);
}



/** \file DoubleValueFunctions.hpp
 * Output functions on a double value.
 * \fn double choose(double n, double k)
 * Find the minimum element.  Usage:
 * \code
 * double theChoose = choose(5, 8);
 * \endcode
 * We output to double since this uses the gamma function.
 */
double choose(double n, double k)
{
    return factorial(n)/(factorial(k)*factorial(n-k));
}



#endif
