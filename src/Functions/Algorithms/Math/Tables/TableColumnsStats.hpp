#ifndef MOLBIOLIB_TABLECOLUMNSSTATS_H
#define MOLBIOLIB_TABLECOLUMNSSTATS_H 1

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
#include "src/Functions/Algorithms/Math/Tables/TableColumnStats.hpp"



/** \file TableColumnsStats.hpp
 * Output statistics on columns of a table.
 * \fn template<size_t M, size_t N, typename... T, typename... U> double tableColumnsPearsonR(Table<T...>& theTable1, Table<U...>& theTable2, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max(), long shift = 0)
 * Compute the Pearson's correlation coefficient.
 * \code
 * Table<int, double, string, double, bool> inputTable1, inputTable2;
 * // Populate inputTable1 and inputTable2.
 * double theCoeff = tableColumnsPearsonR<3, 1>(inputTable1, inputTable2, 0, numeric_limits<size_t>::max(), 0);
 * // 0 means start from 0 and max means to the end (defaults).  The final 0
 * // means to shift by 0.
 * //
 * // One may pass the same table twice, if different columns are compared.
 * \endcode
 * One can find the definition of r at
 * http://en.wikipedia.org/wiki/Correlation_and_dependence
 *
 * Shift Table2's entries +/- by shift relative to Table1 by shift.
 * A plus shift shifts theTable2's entries to the right relative to theTable1.
 * For example, if theTable1 = [0, 1, 1, 1, 0, 0, 0] and
 * theTable2 = [0, 0, 0, 1, 1, 1, 0], then a shift of -2 (not +2) will
 * result in a correlation coefficient of +1 whilst a shift of +1 will
 * result in the comparison of vectors [1, 1, 1, 0, 0, 0] and
 * [0, 0, 0, 1, 1, 1] resulting in a correlation coefficient of -1.
 */
template<size_t M, size_t N, typename... T, typename... U>
double tableColumnsPearsonR(Table<T...>& theTable1, Table<U...>& theTable2, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max(), long shift = 0)
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<M, cols_type>::type& theX = get<M>(theTable1.theData);
    typename tuple_element<N, cols_type>::type& theY = get<N>(theTable2.theData);

#ifdef DEBUG
    if (theTable1.size() == 0)
    {
        cerr << "Error in tableColumnsPearson.  The first input table has no elements, so the coefficient is undefined.  Returning 0.0." << endl;
        return 0.0;
    }
    else if (theTable2.size() == 0)
    {
        cerr << "Error in tableColumnsPearson.  The second input table has no elements, so the coefficient is undefined.  Returning 0.0." << endl;
        return 0.0;
    }
#endif

    size_t shiftSize_t = static_cast<size_t>(labs(shift));

#ifdef DEBUG
    // Either shift is specified, or begin is, but not both.
    assert((shiftSize_t == 0) || (theBegin == 0));
    assert(theX.size() > shiftSize_t);
    assert(theY.size() > shiftSize_t);

    // Below: either user specified theEnd (so assumed subset) or should be
    // long enough after shifting;
    assert((theEnd != numeric_limits<size_t>::max()) ||
           ((shift >= 0) && (theX.size() <= (theY.size() + shiftSize_t))) ||
           ((shift < 0) && (theX.size() >= (theY.size() - shiftSize_t))));
#endif

    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theX.size() - 1;
    }

    double ymean = 0.0, ystddev = 0.0;
    bool shiftPositive = (shift >= 0) ? true : false;
    // true below in StdDev means use sample stddev (1/(N-1))
    if (shiftPositive)
    {
        ymean = tableColumnMean<N>(theTable2, theBegin, theEnd - shiftSize_t);
        ystddev = tableColumnStdDev<N>(theTable2, theBegin, theEnd - shiftSize_t, true);
        theBegin += shiftSize_t;
    }
    else
    {
        ymean = tableColumnMean<N>(theTable2, theBegin + shiftSize_t, theEnd);
        ystddev = tableColumnStdDev<N>(theTable2, theBegin + shiftSize_t, theEnd, true);
        theEnd -= shiftSize_t;
    }
    double xmean = tableColumnMean<M>(theTable1, theBegin, theEnd);
    double xstddev = tableColumnStdDev<M>(theTable1, theBegin, theEnd, true);

    double ssxy = 0.0;
    double tempx, tempy;
    for (size_t i = theBegin; i <= theEnd; ++i)
    {
        // Doing it this way since if there are large values, lose precision
        // if do summation first.
        tempx = static_cast<double>(theX[i]) - xmean;
        if (shiftPositive)
        {
            tempy = static_cast<double>(theY[i - shiftSize_t]) - ymean;
        }
        else
        {
            tempy = static_cast<double>(theY[i + shiftSize_t]) - ymean;
        }
        ssxy += tempx*tempy;
    }

#ifdef DEBUG
    if (xstddev == 0.0)
    {
        cerr << "Error in tableColumnsPearsonR2!  xstddev is zero, so result is infinite.  Returning 0.0 instead." << endl;
        return 0.0;
    }
    if (ystddev == 0.0)
    {
        cerr << "Error in tableColumnsPearsonR2!  ystddev is zero, so result is infinite.  Returning 0.0 instead." << endl;
        return 0.0;
    }
#endif

    double n = static_cast<double>(theEnd - theBegin + 1);

    return (ssxy/((n-1.0)*xstddev*ystddev));
}




/** \file TableColumnsStats.hpp
 * Output statistics on columns of a table.
 * \fn template<size_t M, size_t N, typename... T, typename... U> double tableColumnsPearsonR2(Table<T...>& theTable1, Table<U...>& theTable2, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max(), long shift = 0)
 * Compute the square of the Pearson's correlation coefficient.
 */
template<size_t M, size_t N, typename... T, typename... U>
double tableColumnsPearsonR2(Table<T...>& theTable1, Table<U...>& theTable2, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max(), long shift = 0)
{
    double result = tableColumnsPearsonR<M, N>(theTable1, theTable2, theBegin, theEnd, shift);
    return result*result;
}







/** \file TableColumnsStats.hpp
 * Output statistics on columns of a table.
 * \fn template<size_t M, size_t N, typename... T, typename... U> double tableColumnsSpearman(Table<T...>& theTable1, Table<U...>& theTable2, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max(), bool assumeDiff = false, long shift = 0)
 * Compute the Spearman's correlation coefficient.
 * \code
 * Table<int, double, string, double, bool> inputTable1, inputTable2;
 * // Populate inputTable1 and inputTable2.
 * double theCoeff = tableColumnsSpearman<3, 1>(inputTable1, inputTable2, 0, numeric_limits<size_t>::max(), false, 0);
 * // 0 means start from 0 and max() means to the end (defaults).
 * // The false means do not assume the elements are unique.  (True, speeds
 * // things up some.)  The final 0 means shift by 0.
 * // One may pass the same table twice, if different columns are compared.
 * \endcode
 * One can find the definition of r at
 * http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
 *
 * Shift Table2's entries +/- by shift relative to Table1 by shift.
 * A plus shift shifts theTable2's entries to the right relative to theTable1.
 * For example, if theTable1 = [0, 1, 1, 1, 0, 0, 0] and
 * theTable2 = [0, 0, 0, 1, 1, 1, 0], then a shift of -2 (not +2) will
 * result in a correlation coefficient of +1, since then comparing
 * theTable2 = [0, 1, 1, 1, 0,] vs theTable1 = [0, 1, 1, 1, 0, 0, 0] whilst
 * a shift of +1 will result in the comparison of vectors [1, 1, 1, 0, 0, 0]
 * and [0, 0, 0, 1, 1, 1] resulting in a correlation coefficient of -1.
 */
template<size_t M, size_t N, typename... T, typename... U>
double tableColumnsSpearman(Table<T...>& theTable1, Table<U...>& theTable2, size_t theBegin = 0, size_t theEnd = numeric_limits<size_t>::max(), bool assumeDiff = false, long shift = 0)
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<M, cols_type>::type& theX = get<M>(theTable1.theData);
    typename tuple_element<N, cols_type>::type& theY = get<N>(theTable2.theData);

#ifdef DEBUG
    if (theTable1.size() < 2)
    {
        cerr << "Error in tableColumnsSpearman.  The first input table has no or one elements, so the coefficient is undefined.  Returning 0.0." << endl;
        return 0.0;
    }
    else if (theTable2.size() < 2)
    {
        cerr << "Error in tableColumnsSpearman.  The second input table has no or one elements, so the coefficient is undefined.  Returning 0.0." << endl;
        return 0.0;
    }
#endif

    size_t shiftSize_t = static_cast<size_t>(labs(shift));

#ifdef DEBUG
    // Either shift is specified, or begin is, but not both.
    assert((shiftSize_t == 0) || (theBegin == 0));
    assert(theX.size() > shiftSize_t);
    assert(theY.size() > shiftSize_t);

    // Below: either user specified theEnd (so assumed subset) or should be
    // long enough after shifting;
    assert((theEnd != numeric_limits<size_t>::max()) ||
           ((shift >= 0) && (theX.size() <= (theY.size() + shiftSize_t))) ||
           ((shift < 0) && (theX.size() >= (theY.size() - shiftSize_t))));
#endif

    if (theEnd == numeric_limits<size_t>::max())
    {
        theEnd = theX.size() - 1;
    }

    bool shiftPositive = (shift >= 0) ? true : false;

    // Need to copy the variables in here so can sort and rank without
    // destroying original Tables' values.
    //      X       Y     rank X  rank Y
    Table<double, double, double, double> values;
    tuple<double, double, double, double> aRow;
    // Fill in ranks later...
    get<2>(aRow) = 0.0;
    get<3>(aRow) = 0.0;
    vector<double> X, Y;
    bool breakNow = false;
    for (size_t i = theBegin; !breakNow && (i <= theEnd); ++i)
    {
        if (shift == 0)
        {
            get<0>(aRow) = static_cast<double>(theX[i]);
            get<1>(aRow) = static_cast<double>(theY[i]);
            values.push_back(aRow);
        }
        else if (shiftPositive)
        {
            if ((i+shiftSize_t) > theEnd)
            {
                breakNow = true;
            }
            else
            {
                get<0>(aRow) = static_cast<double>(theX[i+shiftSize_t]);
                get<1>(aRow) = static_cast<double>(theY[i]);
                values.push_back(aRow);
            }
        }
        else
        {
            // shiftPositive is false
            if ((i+shiftSize_t) > theEnd)
            {
                breakNow = true;
            }
            else
            {
                get<0>(aRow) = static_cast<double>(theX[i]);
                get<1>(aRow) = static_cast<double>(theY[i+shiftSize_t]);
                values.push_back(aRow);
            }
        }
    }


    size_t n = values.size();


    // true below means we can sort values.
    tableColumnRank<0, 2>(values, values, assumeDiff, true);
    tableColumnRank<1, 3>(values, values, assumeDiff, true);


    // Now, the last two columns of values have the ranks in them.
    double rho = 0.0;
    double doubleN = static_cast<double>(n);
    if (!assumeDiff)
    {
        // temp = sum((x_i - y_i)^2), where x_i and y_i are the ranks.
        double temp = 0.0, temp2 = 0.0;
        for (size_t i = 0; i < n; ++i)
        {
            temp2 = values.getElement<2>(i) - values.getElement<3>(i);
            temp += temp2*temp2;
        }
        rho = 1.0 - (6.0*temp/(doubleN*(doubleN*doubleN-1.0)));
    }
    else
    {
        double xbar = tableColumnMean<2>(values);
        double ybar = tableColumnMean<3>(values);
        double numer = 0.0, denomX = 0.0, denomY = 0.0;
        double tempX = 0.0, tempY = 0.0;
        for (size_t i = 0; i < n; ++i)
        {
            tempX = values.getElement<2>(i) - xbar;
            tempY = values.getElement<3>(i) - ybar;
            numer += tempX*tempY;
            denomX += tempX*tempX;
            denomY += tempY*tempY;
        }
        rho = numer/sqrt(denomX*denomY);
    }


    return rho;
}



#endif


