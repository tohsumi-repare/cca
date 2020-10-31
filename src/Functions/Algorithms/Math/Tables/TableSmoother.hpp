#ifndef MOLBIOLIB_TABLESMOOTHER_H
#define MOLBIOLIB_TABLESMOOTHER_H 1

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


/** \file TableSmoother.hpp
 * Classes to perform smoothing on a Table's column.
 * Various smoothers are implemented here, all based on quadrature.
 */

#include "src/Objects/Table.hpp"


#ifndef SMOOTHING_EPSILON
// Below may not necessarily be numeric_limits<double>::epsilon(), as
// that is perhaps too small.  Note this is multiplied by the integrand,
// which may be much larger than 1, so decrease a few orders of magnitude.
// #define SMOOTHING_EPSILON numeric_limits<double>::epsilon()
#define SMOOTHING_EPSILON 1.0e-20
#endif


// The derived kernels must all be symmetric with repsect to position around 0.
class SymmetricIntegrandKernel
{
public:
    SymmetricIntegrandKernel()
    {
        border = 0;
    }
    size_t getBorder()
    {
        return border;
    }
    Table<double>& getArray()
    {
        return kernelArray;
    }
protected:
    // Below should never be called.
    virtual double kernel(size_t t)
    {
        assert(true == false);
        return 0.0;
    }
    // Stop when reach this many times epsilon;
    void initializeKernelArray(size_t numEps = 1)
    {
        bool stop = false;
        size_t t = 0;
        double temp = 0.0;
        size_t foundEps = 0;
        while (!stop)
        {
            temp = kernel(t);
            kernelArray.push_back(temp);
            if (fabs(temp) < epsilon)
            {
                // Want to see that the values are decreasing in size.
                // Check last two values.
                if ((kernelArray.size() < 2) ||
                        (fabs(kernelArray[kernelArray.size()-1]) <=
                         fabs(kernelArray[kernelArray.size()-2])))
                {
                    ++foundEps;
                }
                if (foundEps >= numEps)
                {
                    stop = true;
                    border = t;
                    // Push one more extra so that in case have to extend due to
                    // odd-sized interval.
                    ++t;
                    temp = kernel(t);
                    kernelArray.push_back(temp);
                }
            }
            ++t;
        }
    }
    double epsilon;
private:
    Table<double> kernelArray;
    size_t border;
};


// Below is 1/sqrt(2*pi)
#ifndef M_1_SQRT2PI
#define M_1_SQRT2PI 0.3989422804014326779399460599343819L
#endif

class GaussianKernel : public SymmetricIntegrandKernel
{
public:
    GaussianKernel(double theSigma, double theEpsilon = SMOOTHING_EPSILON)
    {
        sigma = theSigma;
        sigma2 = sigma*sigma;
        epsilon = theEpsilon;
        initializeKernelArray();
    }
protected:
    double kernel(size_t t)
    {
        double theT2 = static_cast<double>(t*t);
        return (M_1_SQRT2PI/sigma)*exp(-theT2/(2.0*sigma2));
    }
private:
    double sigma;
    double sigma2;   // sigma^2.
};


class GaussianSecondDerivKernel : public SymmetricIntegrandKernel
{
public:
    GaussianSecondDerivKernel(double theSigma, double theEpsilon = SMOOTHING_EPSILON)
    {
        sigma = theSigma;
        sigma2 = sigma*sigma;
        epsilon = theEpsilon;
        initializeKernelArray();
    }
protected:
    double kernel(size_t t)
    {
        double theT2 = static_cast<double>(t*t);
        return M_1_SQRT2PI*((theT2/sigma2 - 1.0)/(sigma2*sigma))*exp(-theT2/(2.0*sigma2));
    }
private:
    double sigma;
    double sigma2;   // sigma^2.
};



// Below is 1/(sqrt(3)*(pi^(1/4)))
#ifndef M_1_SQRT3PI14
#define M_1_SQRT3PI14 0.8673250705840776348516849306987600L
#endif

class WaveletKernel : public SymmetricIntegrandKernel
{
public:
    WaveletKernel(double theA, double theEpsilon = SMOOTHING_EPSILON)
    {
        a = theA;
        epsilon = theEpsilon;
        initializeKernelArray(2);  // Must hit < |epsilon| twice.
    }
protected:
    double kernel(size_t t)
    {
        double pos = static_cast<double>(t)/a;
        double t2 = pos*pos;
        return M_1_SQRT3PI14*(1.0 - t2)*exp(-t2/2.0);
    }
private:
    double a;
};




/** \file TableSmoother.hpp
 * \fn template<size_t N, typename...T> simpsonIntegral(SymmetricIntegrandKernel& K, Table<T...>& theTable, size_t currRow)
 * Smooths values in a table.
 */
template<size_t N, typename...T>
double simpsonIntegral(SymmetricIntegrandKernel& K, Table<T...>& theTable,
                       size_t currRow)
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    size_t numRows = theColumn.size();
    size_t rowStart = 0, rowEnd = currRow + K.getBorder();
    if (currRow > K.getBorder())
    {
        rowStart = currRow - K.getBorder();
    }
    if ((theColumn.size() - 1) < rowEnd)
    {
        rowEnd = theColumn.size() - 1;
    }

    // If the number of points is not even, make it even.
    if ((rowEnd - rowStart) % 2 != 1)
    {
        if (rowStart > 0)
        {
            --rowStart;
        }
        else
        {
            if (rowEnd < (numRows - 1))
            {
                ++rowEnd;
            }
            else
            {
                --rowEnd;
#ifdef DEBUG
                cerr << "Warning in simpsonIntegral.  Had to subtract 1 from row end to make even number of points for composite Simpson.  May lead to loss of accuracy." << endl;
#endif
            }
        }
    }

    Table<double>& kernelArray = K.getArray();

    // currRow always is >= rowStart.
    double kernelValue = kernelArray[currRow - rowStart];
    double result = static_cast<double>(theColumn[rowStart])*kernelValue;

    // currRow always is <= rowEnd
    kernelValue = kernelArray[rowEnd - currRow];
    result += static_cast<double>(theColumn[rowEnd])*kernelValue;

    size_t j = 1;
    for (size_t i = rowStart + 1; i < rowEnd; ++i, ++j)
    {
        kernelValue = (currRow > i) ? kernelArray[currRow - i] : kernelArray[i - currRow];
        if (j % 2 == 0)
        {
            result += 2.0*static_cast<double>(theColumn[i])*kernelValue;
        }
        else
        {
            result += 4.0*static_cast<double>(theColumn[i])*kernelValue;
        }
    }
    result /= 3.0;

    return result;
}


/** \file TableSmoother.hpp
 * \fn template<size_t N, typename...T> void simpsonIntegralTable(SymmetricIntegrandKernel& K, Table<T...>& theTable, Table<double>& result, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Smooth a column of Table<...double...>.
 * Usage:
 * \code
 * Table<string, double, int, int> test;
 * // fill test array
 * GaussianKernel K(5.0);   // sigma = 5.0;
 * Table<double> smoothedTable;
 * simpsonIntegralTable<0>(K, test, smoothedTable);
 *   // or smoothedTable = simpsonIntegralTable(K, test);
 *   // The above fills smoothedTable with the values.
 * \endcode
 * More generally, one may pass any table:
 * \code
 * Table<long, double, double, string> smoothedTable2;
 * // Below is optional.
 * smoothedTable2.resize(test.size());  // smoothedTable2, must be at least as
 *                                      // long as test.
 * simpsonIntegralTable<0, 2>(K, test, smoothedTable2);
 * // The 2 in the template parameter means to replace smoothedTable2's third
 * // column by the output of simpsonIntegralTable.
 * \endcode
 */
template<size_t N, typename...T>
void simpsonIntegralTable(SymmetricIntegrandKernel& K,
                          Table<T...>& theTable,
                          Table<double>& result,
                          size_t theBegin = 0,
                          size_t theEnd = numeric_limits<size_t>::max())
{
    size_t numRows = theTable.size();
    for (size_t i = 0; i < numRows; ++i)
    {
        result.push_back(simpsonIntegral<N>(K, theTable, i));
    }
}


/** \file TableSmoother.hpp
 * \fn template<size_t N, typename...T> Table<double> simpsonIntegralTable(SymmetricIntegrandKernel& K, Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 */
template<size_t N, typename...T>
Table<double> simpsonIntegralTable(SymmetricIntegrandKernel& K,
                                   Table<T...>& theTable,
                                   size_t theBegin = 0,
                                   size_t theEnd = numeric_limits<size_t>::max())
{
    Table<double> result;
    simpsonIntegralTable<N>(K, theTable, result, theBegin, theEnd);
    return result;
}


/** \file TableSmoother.hpp
 * \fn template<size_t N, size_t M, typename...T, typename...U> void simpsonIntegralTable(SymmetricIntegrandKernel& K, Table<T...>& theTable, Table<U...>& result, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 */
template<size_t N, size_t M, typename...T, typename...U>
void simpsonIntegralTable(SymmetricIntegrandKernel& K,
                          Table<T...>& theTable,
                          Table<U...>& result,
                          size_t theBegin = 0,
                          size_t theEnd = numeric_limits<size_t>::max())
{
    size_t numRows = theTable.size();
    typename Table<U...>::row_type aRow;
    for (size_t i = 0; i < numRows; ++i)
    {
        double resultValue = simpsonIntegral<N>(K, theTable, i);
        result.getRow(i, aRow);
        get<M>(aRow) = resultValue;
        result.setRow(i, aRow);
    }
}






/** \file TableSmoother.hpp
 * \fn template<size_t N, typename...T> void movingAverageTableSmoother(Table<T...>& theTable, size_t halfWindow, Table<double>& result, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 */
template<size_t N, typename...T>
void movingAverageTableSmoother(Table<T...>& theTable,
                                size_t halfWindow,
                                Table<double>& result,
                                size_t theBegin = 0,
                                size_t theEnd = numeric_limits<size_t>::max())
{
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theTable.size() - 1;
    }

    COL_TYPE_TABLE_TYPENAME(N, Table<T...>)& theColumn = get<N>(theTable.theData);

    result.resize(theTable.size(), 0.0);
    size_t loopBegin = theBegin;
    if (theBegin < halfWindow)
    {
        for (size_t i = theBegin; i < halfWindow - theBegin; ++i)
        {
            result[i] = static_cast<double>(theColumn[i]);
        }
        loopBegin = halfWindow - theBegin;
    }

    size_t loopEnd = theEnd;
    if (loopEnd >= (theTable.size() - halfWindow))
    {
        loopEnd = theTable.size() - 1 - halfWindow;
    }

    size_t num = 2*halfWindow + 1;
    for (size_t i = loopBegin; i <= loopEnd; ++i)
    {
        double avg = 0.0;
        for (size_t j = i - halfWindow; j <= i + halfWindow; ++j)
        {
            avg += static_cast<double>(theColumn[j]);
        }
        avg /= static_cast<double>(num);
        result[i] = avg;
    }

    if (loopEnd != theEnd)
    {
        for (size_t i = loopEnd + 1; i <= theEnd; ++i)
        {
            result[i] = static_cast<double>(theColumn[i]);
        }
    }
}




/** \file TableSmoother.hpp
 * \fn template<size_t N, typename...T> Table<double> movingAverageTableSmoother(Table<T...>& theTable, size_t halfWindow, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 */
template<size_t N, typename...T>
Table<double> movingAverageTableSmoother(Table<T...>& theTable,
        size_t halfWindow,
        size_t theBegin = 0,
        size_t theEnd = numeric_limits<size_t>::max())
{
    Table<double> result;
    movingAverageTableSmoother<N>(theTable, halfWindow, result, theBegin, theEnd);
    return result;
}



/** \file TableSmoother.hpp
 * \fn template<size_t N, size_t M, typename...T, typename...U> void movingAverageTableSmoother(Table<T...>& theTable, Table<U...>& result, size_t halfWindow, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 */
template<size_t N, size_t M, typename...T, typename...U>
void movingAverageTableSmoother(Table<T...>& theTable,
                                Table<U...>& result,
                                size_t halfWindow,
                                size_t theBegin = 0,
                                size_t theEnd = numeric_limits<size_t>::max())
{
#ifdef DEBUG
    cerr << "Warning.  movingAverageTableSmoother(theTable, theResult, halfWindow, theBegin, theEnd) is not memory-efficient.  Proceeding." << endl;
#endif
    Table<double> tempResult;
    movingAverageTableSmoother<N>(theTable, halfWindow, result, theBegin, theEnd);
    size_t end = min(theTable.size() - 1, theEnd);
    typename Table<U...>::row_type aRow;
    for (size_t i = theBegin; i <= end; ++i)
    {
        result.getRow(i, aRow);
        get<M>(aRow) = tempResult[i];
        result.setRow(i, aRow);
    }
}


#endif

