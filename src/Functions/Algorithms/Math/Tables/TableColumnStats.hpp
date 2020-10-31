#ifndef MOLBIOLIB_TABLECOLUMNSTATS_H
#define MOLBIOLIB_TABLECOLUMNSTATS_H 1

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


#include "src/Objects/Table.hpp"



/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> size_t tableColumnNumZeros(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Find the number of zeros in a column.
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * size_t numZeros = tableColumnNumZeros<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 */
template<size_t N, typename... T>
size_t tableColumnNumZeros(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
    if (theTable.size() == 0)
    {
        return 0;
    }
    size_t numZeros = 0;
    for (size_t i = theBegin; i <= theEnd; ++i)
    {
        if (static_cast<double>(theColumn[i]) == 0.0)
        {
            ++numZeros;
        }
    }

    return numZeros;
}


/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnMin(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Find the minimum element.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double theMin = tableColumnMin<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 */
template<size_t N, typename... T>
double tableColumnMin(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
#ifdef DEBUG
    if (theTable.size() == 0)
    {
        cerr << "Error in tableColumnMin.  The input table has no elements, so the minimum is undefined.  Returning 0.0." << endl;
        return 0.0;
    }
#endif
    double theMin = static_cast<double>(theColumn[theBegin]);
    for (size_t i = theBegin+1; i <= theEnd; ++i)
    {
        if (static_cast<double>(theColumn[i]) < theMin)
        {
            theMin = static_cast<double>(theColumn[i]);
        }
    }

    return theMin;
}



/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnMax(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Find the maximum element.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double theMax = tableColumnMax<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 */
template<size_t N, typename... T>
double tableColumnMax(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
#ifdef DEBUG
    if (theTable.size() == 0)
    {
        cerr << "Error in tableColumnMax.  The input table has no elements, so the minimum is undefined.  Returning 0.0." << endl;
        return 0.0;
    }
#endif
    double theMax = static_cast<double>(theColumn[theBegin]);
    for (size_t i = theBegin+1; i <= theEnd; ++i)
    {
        if (theMax < static_cast<double>(theColumn[i]))
        {
            theMax = static_cast<double>(theColumn[i]);
        }
    }

    return theMax;
}



/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnSum(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Compute the sum.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double theSum = tableColumnSum<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 */
template<size_t N, typename... T>
double tableColumnSum(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
    if (theColumn.size() == 0)
    {
#ifdef DEBUG
        cerr << "Warning in tableColumnSum.  The input table has no elements, so the sum is undefined.  Returning 0.0." << endl;
#endif
        return 0.0;
    }
    double theSum = 0.0;
    for (size_t i = theBegin; i <= theEnd; ++i)
    {
        theSum += static_cast<double>(theColumn[i]);
    }

    return theSum;
}



/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnMean(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Compute the mean (average).  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double theMean = tableColumnMean<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 */
template<size_t N, typename... T>
double tableColumnMean(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
    if (theColumn.size() == 0)
    {
#ifdef DEBUG
        cerr << "Warning in tableColumnMean.  The input table has no elements, so the mean is undefined.  Returning 0.0." << endl;
#endif
        return 0.0;
    }
    double theMean = tableColumnSum<N>(theTable, theBegin, theEnd);
    theMean = theMean/static_cast<double>(theEnd - theBegin + 1);

    return theMean;
}



/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnStdDev(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max(), bool isSample = false)
 * Compute the standard deviation.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double theStdDev = tableColumnStdDev<2>(inputTable, 0, numeric_limits<size_t>::max(), false);
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * // false means do not use sample stdDev which normalizes with 1/(N-1)
 * // instead of 1/N.
 * \endcode
 */
template<size_t N, typename... T>
double tableColumnStdDev(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max(), bool isSample = false)
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
    if (theColumn.size() == 0)
    {
#ifdef DEBUG
        cerr << "Warning in tableColumnStdDev.  The input table has no elements, so the standard deviation is undefined.  Returning 0.0." << endl;
#endif
        return 0.0;
    }
    double theMean = tableColumnMean<N>(theTable, theBegin, theEnd);
    double sqrStdDev = 0.0;
    for (size_t i = theBegin; i <= theEnd; ++i)
    {
        sqrStdDev += (static_cast<double>(theColumn[i]) - theMean)*
                     (static_cast<double>(theColumn[i]) - theMean);
    }
    // If sample then, normalize 1/(N-1) else normalize 1/N.
    // Since equality in theEnd, should have +1 below if normalize 1/N.
    size_t sampleAdjust = isSample ? 0 : 1;
    sqrStdDev = sqrStdDev/static_cast<double>(theEnd - theBegin + sampleAdjust);

    return sqrt(sqrStdDev);
}


/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnMedian(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Compute the median.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double theMedian = tableColumnMedian<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 */
template<size_t N, typename... T>
double tableColumnMedian(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);

    Table<double> tempVector;
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
    if (theColumn.size() == 0)
    {
#ifdef DEBUG
        cerr << "Warning in tableColumnMedian.  The input table has no elements, so the mean is undefined.  Returning 0.0." << endl;
#endif
        return 0.0;
    }
    for (size_t i = theBegin; i <= theEnd; ++i)
    {
        tempVector.push_back(static_cast<double>(theColumn[i]));
    }

    sortTable<0>(tempVector);
    size_t vecSize = tempVector.size();

    if ((vecSize % 2) == 1)
    {
        return tempVector[(vecSize - 1)/2];
    }
    else
    {
        size_t tempIndex = (vecSize + 1)/2;
        return (tempVector[tempIndex] + tempVector[tempIndex-1])/2.0;
    }
}


/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> tuple<double, double, double> tableColumnQ123(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Compute the mean (average).  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * tuple<double, double, double> q1q2q3 = tableColumnQ123<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 * Based on
 * https://stackoverflow.com/questions/11964552/finding-quartiles
 * Using the Matlab version.
 */
template<size_t N, typename... T>
tuple<double, double, double> tableColumnQ123(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);

    tuple<double, double, double> result(0.0, 0.0, 0.0);
    Table<double> tempVector;
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
    if (theColumn.size() == 0)
    {
#ifdef DEBUG
        cerr << "Warning in tableColumnMedian.  The input table has no elements, so the mean is undefined.  Returning 0.0." << endl;
#endif
        return result;
    }
    for (size_t i = theBegin; i <= theEnd; ++i)
    {
        tempVector.push_back(static_cast<double>(theColumn[i]));
    }

    sortTable<0>(tempVector);
    size_t vecSize = tempVector.size();
    if (vecSize == 1) {
      get<0>(result) = tempVector[0];
      get<1>(result) = tempVector[0];
      get<2>(result) = tempVector[0];
    } else {
       // lerp(x, y, z) = (1 - z)*x + z*y;
       // firstQuantile
       // double poi = lerp(-0.5, vecSize - 0.5, 0.25);
       double poi = (1.0 - 0.25)*(-0.5) + 0.25*(vecSize - 0.5);
       size_t left = max(static_cast<size_t>(floor(poi)),
                         static_cast<size_t>(0));
       size_t right = min(static_cast<size_t>(ceil(poi)), vecSize - 1);
       double dataLeft = tempVector[left];
       double dataRight = tempVector[right];
       // get<0>(result) = lerp(dataLeft, dataRight, poi - left);
       get<0>(result) = (1.0 - (poi - left))*dataLeft + (poi - left)*dataRight;
       // median
       // poi = lerp(-0.5, vecSize - 0.5, 0.5);
       poi = (1.0 - 0.5)*(-0.5) + 0.5*(vecSize - 0.5);
       left = max(static_cast<size_t>(floor(poi)), static_cast<size_t>(0));
       right = min(static_cast<size_t>(ceil(poi)), vecSize - 1);
       dataLeft = tempVector[left];
       dataRight = tempVector[right];
       // get<1>(result) = lerp(dataLeft, dataRight, poi - left);
       get<1>(result) = (1.0 - (poi - left))*dataLeft + (poi - left)*dataRight;
       // thirdQuantile
       // poi = lerp(-0.5, vecSize - 0.5, 0.75);
       poi = (1.0 - 0.75)*(-0.5) + 0.75*(vecSize - 0.5);
       left = max(static_cast<size_t>(floor(poi)), static_cast<size_t>(0));
       right = min(static_cast<size_t>(ceil(poi)), vecSize - 1);
       dataLeft = tempVector[left];
       dataRight = tempVector[right];
       // get<2>(result) = lerp(dataLeft, dataRight, poi - left);
       get<2>(result) = (1.0 - (poi - left))*dataLeft + (poi - left)*dataRight;
    }
    return result;
}


/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnStdDevFromMedian(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Compute the standard deviation from the median.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double theStdDev = tableColumnStdDevFromMedian<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 */
template<size_t N, typename... T>
double tableColumnStdDevFromMedian(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
    if (theColumn.size() == 0)
    {
#ifdef DEBUG
        cerr << "Warning in tableColumnStdDevFromMedian.  The input table has no elements, so the standard deviation is undefined.  Returning 0.0." << endl;
#endif
        return 0.0;
    }
    double theMedian = tableColumnMedian<N>(theTable, theBegin, theEnd);
    double sqrStdDev = 0.0;
    for (size_t i = theBegin; i <= theEnd; ++i)
    {
        sqrStdDev += (static_cast<double>(theColumn[i]) - theMedian)*
                     (static_cast<double>(theColumn[i]) - theMedian);
    }
    sqrStdDev = sqrStdDev/static_cast<double>(theEnd - theBegin + 1);

    return sqrt(sqrStdDev);
}


/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnMeanAbsDev(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Compute the mean (average).  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double meanAbsDev = tableColumnMeanAbsDev<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 */
template<size_t N, typename... T>
double tableColumnMeanAbsDev(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);

    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
    if (theColumn.size() == 0)
    {
#ifdef DEBUG
        cerr << "Warning in tableColumnAbsDev.  The input table has no elements, so the deviation is undefined.  Returning 0.0." << endl;
#endif
        return 0.0;
    }

   double median = tableColumnMedian<N>(theTable, theBegin, theEnd);
   Table<double> result;

   for (size_t i = theBegin; i <= theEnd; ++i)
   {
       result.push_back(fabs(theColumn[i] - median));
   }

   return tableColumnMedian<0>(result);
}



/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnSkew(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Compute the skew.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double theSkew = tableColumnSkew<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 */
template<size_t N, typename... T>
double tableColumnSkew(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
    if (theColumn.size() == 0)
    {
#ifdef DEBUG
        cerr << "Warning in tableColumnSkew.  The input table has no elements, so the skew is undefined.  Returning 0.0." << endl;
#endif
        return 0.0;
    }

    double theMean = tableColumnMean<N>(theTable, theBegin, theEnd);

    double mu3 = 0.0, sigma2 = 0.0;
    for (size_t i = theBegin; i <= theEnd; ++i)
    {
        double temp = (theColumn[i] - theMean);
        double t2 = temp*temp;
        mu3 += temp*t2;
        sigma2 += t2;
    }
    double n = static_cast<double>(theEnd - theBegin + 1);
    mu3 = mu3/n;
    sigma2 = sigma2/n;
    sigma2 = sqrt(sigma2);

    return mu3/(sigma2*sigma2*sigma2);
}



/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnKurtosis(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Compute the kurtosis.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double theKurtosis = tableColumnKurtosis<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 */
template<size_t N, typename... T>
double tableColumnKurtosis(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }
    if (theColumn.size() == 0)
    {
#ifdef DEBUG
        cerr << "Warning in tableColumnKurtosis.  The input table has no elements, so the kurtosis is undefined.  Returning 0.0." << endl;
#endif
        return 0.0;
    }

    double theMean = tableColumnMean<N>(theTable, theBegin, theEnd);

    double mu4 = 0.0, sigma2 = 0.0;
    for (size_t i = theBegin; i <= theEnd; ++i)
    {
        double temp = (theColumn[i] - theMean);
        double t2 = temp*temp;
        mu4 += t2*t2;
        sigma2 += t2;
    }
    double n = static_cast<double>(theEnd - theBegin + 1);
    mu4 = mu4/n;
    sigma2 = sigma2/n;

    return (mu4/(sigma2*sigma2)) - 3.0;
}



/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnJarqueBera(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Compute the Jarque-Bera test for normality (0 = normal).  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double theJarqueBera = tableColumnJarqueBera<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 * Note that the Jarque-Bera statistic gets more sensitive as the number of
 * entries, n, gets larger.
 */
template<size_t N, typename... T>
double tableColumnJarqueBera(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{

    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }

    size_t numRows = theEnd - theBegin + 1;
#ifdef DEBUG
    assert(numRows >= 4);
#endif

    double S = tableColumnSkew<N>(theTable, theBegin, theEnd);
    double K = tableColumnKurtosis<N>(theTable, theBegin, theEnd);
    double n = static_cast<double>(numRows);

    return ((n/6.0)*(S*S + K*K/4.0));
}


/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t N, typename... T> double tableColumnDAgostinoK2(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
 * Compute the D'Agostino K^2 value.  See
 * http://en.wikipedia.org/wiki/D'Agostino's_K-squared_test and
 * D'Agostino, Ralph B., Albert Belanger, and Ralph B. D'Agostino, Jr.
 * "A Suggestion for Using Powerful and Informative Tests of Normality",
 * The American Statistician, Vol. 44, No. 4. (Nov., 1990), pp. 316–321.
 * Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * double theDAgostinoK2 = tableColumnDAgostinoK2<2>(inputTable, 0, numeric_limits<size_t>::max());
 * // 0 means start from 0 and numeric_limits<size_t>::max() means to the end (defaults).
 * \endcode
 * If the null hypothesis of normality is true, the K^2 should be
 * Chi^2 distributed with 2 degrees of freedom.  Now, the Chi cumulative
 * distribution function is gamma(k/2, x/2)/Gamma(k/2), where gamma(s, x) is
 * the incomplete gamma function and Gamma(x) is the gamma function.
 * Since k = 2, and gamma(1, x) = 1 - exp(-x) and Gamma(1/2) = sqrt(pi), we
 * have that the p-value = 1.0 - exp(-1.0*theDAgostinoK2/2.0).
 */
template<size_t N, typename... T>
double tableColumnDAgostinoK2(Table<T...>& theTable, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max())
{

    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theColumn.size() - 1;
    }

    size_t numRows = theEnd - theBegin + 1;
#ifdef DEBUG
    assert(numRows >= 9);
#endif

    double g1 = tableColumnSkew<N>(theTable, theBegin, theEnd);
    double g2 = tableColumnKurtosis<N>(theTable, theBegin, theEnd);
    double n = static_cast<double>(numRows);

    double mu2 = 6.0*(n-2.0)/((n+1.0)*(n+3.0));
    double gamma2 = 36.0*(n-7.0)*(n*n + 2.0*n - 5.0)/
                    ((n-2.0)*(n+5.0)*(n+7.0)*(n+9.0));
    double W2 = sqrt(2.0*gamma2 + 4.0) - 1.0;
    double W = sqrt(W2);
    double delta = 1.0/sqrt(log(W));
    double alpha2 = 2.0/(W2 - 1.0);
    double alpha = sqrt(alpha2);

    double Z1 = delta*log(g1/(alpha*sqrt(mu2)) + sqrt(g1*g1/(alpha2*mu2) + 1.0));


    double mu1 = -6.0/(n+ 1.0);
    mu2 = 24.0*n*(n-2.0)*(n-3.0)/
          ((n+1.0)*(n+1.0)*(n+3.0)*(n+5.0));
    double gamma1 = (6.0*(n*n - 5.0*n + 2.0)/
                     ((n+7.0)*(n+9.0)))
                    *sqrt(6.0*(n+3.0)*(n+5.0)/(n*(n-2.0)*(n-3.0)));
    double A = 6.0 + (8.0/gamma1)*(2.0/gamma1 + sqrt(1.0 + 4.0/(gamma1*gamma1)));

    double Z2 = sqrt(9.0*A/2.0)*(1.0 - 2/(9.0*A) -
                                 pow((1.0-2.0/A)/(1.0+((g2-mu1)/sqrt(mu2))*sqrt(2.0/(A - 4.0))),
                                     1.0/3.0));

    return (Z1*Z1 + Z2*Z2);

}


/** \file TableColumnStats.hpp
 * Output statistics on a column of a table.
 * \fn template<size_t M, size_t N, typename... T, typename... U> void tableColumnRank(Table<T...>& valTable, Table<U...>& rankTable, bool assumeDiff = false, bool canSort = false)
 * Ranks the values in valTable's column M and outputs the ranks in rankTable's
 * column N (both M and N are 0-indexed).  ranks start from 1.  rankTable's
 * column N must be of some numeric type (warning: if not assumeDiff, this
 * should be a floating point type).  If there are not enough elements,
 * rankTable will be resize to valTable's length.  Assuming M != N, rankTable
 * can be the same as valTable.  If assumeDiff is true, then one assumes all
 * values in valTable are unique and thus need not have to average rank.  Even
 * if the values in valTable are not unique and one wants to not average rank
 * same values, then pass assumeDiff = true.  If canSort is true, then sorting
 * will take place in valTable, thus saving memory.  Otherwise, valTable's
 * values will be copied to a temporary table, sorted, and mapped (which is
 * very memory intensive).
 */
template<size_t M, size_t N, typename... T, typename... U>
void tableColumnRank(Table<T...>& valTable, Table<U...>& rankTable,
                     bool assumeDiff = false, bool canSort = false)
{
    typedef typename Table<T...>::row_type row_typeV;
    typedef typename tuple_element<M, row_typeV>::type valType;
    typedef typename Table<U...>::row_type row_typeR;
    typedef typename tuple_element<N, row_typeR>::type rankType;

    size_t n = valTable.size();
    if (rankTable.size() < n)
    {
        rankTable.resize(n);
    }


    if (canSort)
    {
        sortTable<M>(valTable);
        for (size_t i = 0; i < n; ++i)
        {
            // Normally, below would work, except compiler gets confused on
            // if <N is to be parsed as operator<.  So, we do it manually.
            // This is happening because valType and rankType are templated
            // types.
            //     We will denote other places below as "operator< problem".
            // rankTable.setElement<N>(i, static_cast<rankType>(i+1));
            (get<N>(rankTable.theData))[i] = static_cast<rankType>(i+1);
        }
        if (!assumeDiff)
        {
            bool inSame = false;
            size_t i = 0, start = 0, end = 0;
            while (i < (n-1))
            {
                // operator< problem.
                // if (valTable.getElement<M>(i) == valTable.getElement<M>(i+1))
                if ((get<M>(valTable.theData))[i] ==
                        (get<M>(valTable.theData))[i+1])
                {
                    if (!inSame)
                    {
                        start = i;
                        inSame = true;
                    }
                    end = i+1;
                    // Below if is to handle the case of similar values
                    // at the end.
                    if (inSame && (i == (n-2)))
                    {
                        double theMean = tableColumnMean<N>(rankTable,
                                                            start, end);
                        for (size_t j = start; j <= end; ++j)
                        {
                            // operator< problem.
                            // rankTable.setElement<N>(j, static_cast<rankType>(theMean));
                            (get<N>(rankTable.theData))[j] = static_cast<rankType>(theMean);
                        }
                    }
                }
                else
                {
                    if (inSame)
                    {
                        inSame = false;
                        double theMean = tableColumnMean<N>(rankTable,
                                                            start, end);
                        for (size_t j = start; j <= end; ++j)
                        {
                            // operator< problem.
                            // rankTable.setElement<N>(j, static_cast<rankType>(theMean));
                            (get<N>(rankTable.theData))[j] = static_cast<rankType>(theMean);
                        }
                    }
                }
                ++i;
            }
        }
    }
    else
    {
        // Cannot sort, so must copy values to some other table.
        Table<size_t, valType, rankType> values;
        tuple<size_t, valType, rankType> aRow;
        for (size_t i = 0; i < n; ++i)
        {
            get<0>(aRow) = i;
            // operator< problem.
            // get<1>(aRow) = valTable.getElement<M>(i);
            get<1>(aRow) = (get<M>(valTable.theData))[i];
            get<2>(aRow) = 0;    // Fill in rank later.
            values.push_back(aRow);
        }
        sortTable<1>(values);
        for (size_t i = 0; i < n; ++i)
        {
            // operator< problem.
            // values.setElement<2>(i, static_cast<rankType>(i+1));
            (get<2>(values.theData))[i] = static_cast<rankType>(i+1);
        }
        if (!assumeDiff)
        {
            bool inSame = false;
            size_t i = 0, start = 0, end = 0;
            while (i < (n-1))
            {
                // operator< problem.
                // if (values.getElement<1>(i) == values.getElement<1>(i+1))
                if ((get<1>(values.theData))[i] ==
                        (get<1>(values.theData))[i+1])
                {
                    if (!inSame)
                    {
                        start = i;
                        inSame = true;
                    }
                    end = i+1;
                    // Below if is to handle the case of similar values
                    // at the end.
                    if (inSame && (i == (n-2)))
                    {
                        double theMean = tableColumnMean<2>(values,
                                                            start, end);
                        for (size_t j = start; j <= end; ++j)
                        {
                            // operator< problem.
                            // values.setElement<2>(j, static_cast<rankType>(theMean));
                            (get<2>(values.theData))[j] = static_cast<rankType>(theMean);
                        }
                    }
                }
                else
                {
                    if (inSame)
                    {
                        inSame = false;
                        double theMean = tableColumnMean<2>(values,
                                                            start, end);
                        for (size_t j = start; j <= end; ++j)
                        {
                            // operator< problem.
                            // values.setElement<2>(j, static_cast<rankType>(theMean));
                            (get<2>(values.theData))[j] = static_cast<rankType>(theMean);
                        }
                    }
                }
                ++i;
            }
        }

        for (size_t i = 0; i < n; ++i)
        {
            values.getRow(i, aRow);
            // operator< problem.
            // rankTable.setElement<N>(get<0>(aRow), get<2>(aRow));
            (get<N>(rankTable.theData))[get<0>(aRow)] = get<2>(aRow);
        }
    }
}


#endif


