#ifndef MOLBIOLIB_TABLEHISTOGRAM_H
#define MOLBIOLIB_TABLEHISTOGRAM_H 1

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


/** \file TableHistogram.hpp
 * Output a histogram table.
 * \fn template<size_t N, typename... T> void tableHistogram(Table<T...>& theTable, Table<typename tuple_element<N, typename Table<T...>::row_type>::type, size_t>& histogramTable)
 * Generate a histogram.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * // Below example is if we are doing a histogram on the double column:
 * Table<double, size_t> histoTable;
 * tableHistogram<2>(inputTable, histoTable);
 * \endcode
 */
template<size_t N, typename... T>
void tableHistogram(Table<T...>& theTable,
                    Table<typename tuple_element<N, typename Table<T...>::row_type>::type,
                    size_t>& histogramTable)
{
    typedef typename tuple_element<N, typename Table<T...>::row_type>::type HistogramElementType;
    map<HistogramElementType, size_t> theHistogram;
    typename Table<T...>::row_type currentRow;
    HistogramElementType currentElement;
    size_t numRows = theTable.size();
    for (size_t i = 0; i < numRows; ++i)
    {
        theTable.getRow(i, currentRow);
        currentElement = get<N>(currentRow);
        if (theHistogram.find(currentElement) == theHistogram.end())
        {
            theHistogram[currentElement] = 1;
        }
        else
        {
            theHistogram[currentElement] += 1;
        }
    }

    tuple<HistogramElementType, size_t> histogramRow;
    for (typename map<HistogramElementType, size_t>::iterator i = theHistogram.begin();
            i != theHistogram.end(); ++i)
    {
        get<0>(histogramRow) = i->first;
        get<1>(histogramRow) = i->second;
        histogramTable.push_back(histogramRow);
    }

}



/** \file TableHistogram.hpp
 * Output a histogram table.
 * \fn template<size_t N, typename... T> Table<typename tuple_element<N, typename Table<T...>::row_type>::type, size_t> tableHistogram(Table<T...>& theTable)
 * Generate a histogram.  Given the declarations for the pass-by-reference
 * version of #tableHistogram, then the usage of this function is:
 * \code
 * histoTable = tableHistogram<2>(inputTable);
 * \endcode
 */
template<size_t N, typename... T>
Table<typename tuple_element<N, typename Table<T...>::row_type>::type, size_t> tableHistogram(Table<T...>& theTable)
{
    Table<typename tuple_element<N, typename Table<T...>::row_type>::type, size_t> histogramTable;
    tableHistogram<N>(theTable, histogramTable);
    return histogramTable;
}





/** \file TableHistogram.hpp
 * Output a histogram table.
 * \fn template<size_t N, typename... T> void tableDoublesHistogram(Table<T...>& theTable, Table<double, size_t>& histogramTable, size_t numDigits = numeric_limits<size_t>::max(), string rounding = "round", bool inverseData = false)
 * Generate a histogram.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * // Below example is if we are doing a histogram on the double column:
 * Table<double, size_t> histoTable;
 * tableDoublesHistogram<2>(inputTable, histoTable);
 * \endcode
 * Use <code>numDigits = 0</code> if want to round to the nearest integer.
 */
template<size_t N, typename... T>
void tableDoublesHistogram(Table<T...>& theTable,
                           Table<double, size_t>& histogramTable,
                           size_t numDigits = numeric_limits<size_t>::max(),
                           string rounding = "round",
                           bool inverseData = false)
{
    bool doFloor = false, doRound = false, doCeiling = false;
    if (rounding == "round")
    {
        doRound = true;
    }
    else if (rounding == "floor")
    {
        doFloor = true;
    }
    else if (rounding == "ceiling")
    {
        doCeiling = true;
    }
#ifdef DEBUG
    else
    {
        assert(true == false);    // this should never happen.
    }
#endif
    double tenPow = 1.0;
    if (numDigits != numeric_limits<size_t>::max())
    {
        tenPow = pow(10.0, static_cast<double>(numDigits));
    }
    map<double, size_t> theHistogram;
    typename Table<T...>::row_type currentRow;
    double currentElement;
    size_t numRows = theTable.size();
    for (size_t i = 0; i < numRows; ++i)
    {
        theTable.getRow(i, currentRow);
        currentElement = get<N>(currentRow);
        if (inverseData)
        {
            currentElement = 1.0/currentElement;
        }
        if (numDigits != numeric_limits<size_t>::max())
        {
            currentElement *= tenPow;
            if (doRound)
            {
                currentElement = round(currentElement);
            }
            else if (doFloor)
            {
                currentElement = floor(currentElement);
            }
            else if (doCeiling)
            {
                currentElement = ceil(currentElement);
            }
#ifdef PROG_DEBUG
            else
            {
                assert(true == false);    // This should never happen.
            }
#endif
            currentElement /= tenPow;
        }

        if (theHistogram.find(currentElement) == theHistogram.end())
        {
            theHistogram[currentElement] = 1;
        }
        else
        {
            theHistogram[currentElement] += 1;
        }
    }

    tuple<double, size_t> histogramRow;
    for (typename map<double, size_t>::iterator i = theHistogram.begin();
            i != theHistogram.end(); ++i)
    {
        get<0>(histogramRow) = i->first;
        get<1>(histogramRow) = i->second;
        histogramTable.push_back(histogramRow);
    }
}



/** \file TableHistogram.hpp
 * Output a histogram table.
 * \fn template<size_t N, typename... T> Table<double, size_t> tableDoublesHistogram(Table<T...>& theTable, size_t numDigits = numeric_limits<size_t>::max(), string rounding = "round", bool inverseData = false)
 * Generate a histogram.  Given the declarations for the pass-by-reference
 * version of #tableHistogram, then the usage of this function is:
 * \code
 * histoTable = tableDoublesHistogram<2>(inputTable);
 * \endcode
 */
template<size_t N, typename... T>
Table<double, size_t> tableDoublesHistogram(Table<T...>& theTable,
        size_t numDigits = numeric_limits<size_t>::max(),
        string rounding = "round",
        bool inverseData = false)
{
    Table<double, size_t> histogramTable;
    tableDoublesHistogram<N>(theTable, histogramTable, numDigits, rounding, inverseData);
    return histogramTable;
}






#endif


