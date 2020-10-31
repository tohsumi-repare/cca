#ifndef MOLBIOLIB_HYPERGEOMETRIC_H
#define MOLBIOLIB_HYPERGEOMETRIC_H 1

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
#include "src/Functions/Algorithms/Math/Scalar/AppliedStatisticsAlgorithms.hpp"


/** \file Hypergeometric.hpp
 * Hypergeometric distribution
 * \fn double hypergeometricPLessThanOrEqualSampleSuccess(int numPop, int popSuccess, int numSample, int sampleSuccess)
 * This gives P(X <= sampleSuccess).  All input values are <code>int</code>'s
 * because of using the ASA152 routines.
 * \code
 * double theValue = hypergeometricPLessThanOrEqualSampleSuccess(5000, 2000, 50, 20);
 * \endcode
 */
double hypergeometricPLessThanOrEqualSampleSuccess(int numPop,
        int popSuccess,
        int numSample,
        int sampleSuccess)
{
    // Use ASA152 routine.
    int ifault;
    // false = do not want point, but rather CDF.
    double result = chyper(false, numSample, sampleSuccess, numPop, popSuccess, &ifault);
#ifdef DEBUG
    if (ifault != 0)
    {
        cerr << "ifault from chyper of hypergeometricPLessThanOrEqualSampleSuccess = " << ifault << endl;
    }
#endif
    return result;
}

#endif


