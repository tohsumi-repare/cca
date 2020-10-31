#ifndef MOLBIOLIB_POISSONFINDMINNUM_H
#define MOLBIOLIB_POISSONFINDMINNUM_H 1

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


/** \file PoissonFindMinNum.hpp
 * Compute the minimum number to satify a Poisson distribution FDR.
 * \fn size_t poissonFindMinNum(double lambda, double FDR)
 * Finds the first integer value number satisfying the FDR.
 */
size_t poissonFindMinNum(double lambda, double FDR)
{
#ifdef DEBUG
    assert(FDR > 0.0);
#endif
    if (lambda == 0.0)
    {
        return 1;    // No coverage on this contig,
    }
    // so return a MIN_NUM that can't happen.
    size_t result = static_cast<size_t>(floor(lambda));
    double cdf = 1.0;
    int ifault = 0;
    while (cdf > FDR)
    {
        ++result;
        double x = static_cast<double>(result);
        cdf = gammad(lambda, x, &ifault);
    }
    return result;
}


/** \file PoissonFindMinNum.hpp
 * Compute the minimum number to satify a Poisson distribution FDR.
 * \fn size_t poissonFindMinNumDouble(double lambda, double FDR, double precision = 1.0e-8, size_t maxIter = 1000)
 * Finds the approximate closest double value number satisfying the FDR.
 */
double poissonFindMinNumDouble(double lambda, double FDR,
                               double precision = 1.0e-8,
                               size_t maxIter = 1000)
{
    if (lambda == 0.0)
    {
        return 1.0;    // No coverage on this contig,
    }
    // so return a MIN_NUM that can't happen.

    double start = lambda, end = lambda*2.0;
    size_t counter = 0;
    int ifault = 0;

    // First find interval of start and end by doubling inteveral.
    while ((counter < maxIter) &&
            gammad(lambda, end, &ifault) > FDR)
    {
        end *= 2.0;
        ++counter;
    }
#ifdef DEBUG
    if (counter >= maxIter)
    {
        cerr << "Warning in poissonFindMinNumDouble.  Maximum number of iterations reached.  Accuracy may be compromised.  Continuing." << endl;
    }
#endif
    start = end/2.0;


    // Now use bisection method.
    counter = 0;
    while ((counter < maxIter) && (fabs(end - start) > precision))
    {
        double midPoint = (start + end)/2.0;
        double value = gammad(lambda, midPoint, &ifault);
        if (value < FDR)
        {
            end = midPoint;
        }
        else if (value > FDR)
        {
            start = midPoint;
        }
        else
        {
            return midPoint;    // On the rare chance equality is achieved.
        }
        ++counter;
    }

    return (start + end)/2.0;
}



#endif

