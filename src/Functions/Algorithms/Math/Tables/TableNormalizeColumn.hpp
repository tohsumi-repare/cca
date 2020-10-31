#ifndef MOLBIOLIB_TABLENORMALIZECOLUMN_H
#define MOLBIOLIB_TABLENORMALIZECOLUMN_H 1

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


/** \file TableNormalizeColumn.hpp
 * Normalize a column of a table so it's magnitude is 1.
 * \fn template<size_t N, typename... T> void tableEuclideanNormalizeColumn(Table<T...>& theTable)
 * Normalize a column.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * tableEuclideanNormalizeColumn<2>(inputTable);
 * \endcode
 */
template<size_t N, typename... T>
void tableEuclideanNormalizeColumn(Table<T...>& theTable)
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
#ifdef DEBUG
    if (theTable.size() == 0)
    {
        cerr << "Error in tableEuclideanNoramlizeColumn!  The input table has no elements.  Returning the Table unchanged." << endl;
        return;
    }
#endif

    size_t numRows = theTable.size();

    double temp, xbar = 0.0;
    for (size_t i = 0; i < numRows; ++i)
    {
        temp = static_cast<double>(theColumn[i]);
        xbar += temp*temp;
    }
    xbar = sqrt(xbar);
#ifdef DEBUG
    if (xbar == 0.0)
    {
        cerr << "Error in tableEuclideanNoramlizeColumn!  The norm is zero.  Returning the Table unchanged." << endl;
        return;
    }
#endif

    for (size_t i = 0; i < numRows; ++i)
    {
        theColumn[i] /= xbar;
    }
}




/** \file TableNormalizeColumn.hpp
 * Normalize a column of a table so it's absolute valued sum is 1.
 * \fn template<size_t N, typename... T> void table1NormNormalizeColumn(Table<T...>& theTable)
 * Normalize a column.  Usage:
 * \code
 * Table<int, string, double, bool> inputTable;
 * // Populate inputTable
 * table1NormNormalizeColumn<2>(inputTable);
 * \endcode
 */
template<size_t N, typename... T>
void table1NormNormalizeColumn(Table<T...>& theTable)
{
    typedef typename Table<T...>::cols_type cols_type;
    typename tuple_element<N, cols_type>::type& theColumn = get<N>(theTable.theData);
#ifdef DEBUG
    if (theTable.size() == 0)
    {
        cerr << "Error in table1NormNoramlizeColumn!  The input table has no elements.  Returning the Table unchanged." << endl;
        return;
    }
#endif

    size_t numRows = theTable.size();

    double tempSum = 0.0;
    for (size_t i = 0; i < numRows; ++i)
    {
        tempSum += abs(static_cast<double>(theColumn[i]));
    }
#ifdef DEBUG
    if (tempSum == 0.0)
    {
        cerr << "Error in table1NormNoramlizeColumn!  The norm is zero.  Returning the Table unchanged." << endl;
        return;
    }
#endif

    for (size_t i = 0; i < numRows; ++i)
    {
        theColumn[i] /= tempSum;
    }
}




#endif


