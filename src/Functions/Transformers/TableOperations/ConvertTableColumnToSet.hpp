#ifndef MOLBIOLIB_CONVERTTABLECOLUMNTOSET_H
#define MOLBIOLIB_CONVERTTABLECOLUMNTOSET_H 1

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



/** \file ConvertTableColumnToSet.hpp
 * Converts a Table's column into a set.
 * \fn template<size_t N, typename... T> void convertTableColumnToSet(Table<T...>& theTable, set<typename tuple_element<N, Table<T...>::row_type>::type>& outputSet)
 * \section Usage Usage:
 * \code
 * Table<string, Sequence, unsigned long> myExtendedFastaTable;
 * // Add data to myExtendedFastaTable
 * set<Sequence> theSeqs;
 * convertTableColumnToSet<1>(myExtendedFastaTable, theSeqs);
 * \endcode
 */
template<size_t N, typename... T>
void convertTableColumnToSet(Table<T...>& theTable,
                             set<typename tuple_element<N, typename Table<T...>::row_type>::type>& outputSet)
{
    outputSet.clear();
    COL_TYPE_TABLE_TYPENAME(N, Table<T...>)& theColumn = get<N>(theTable.theData);
    size_t numRows = theColumn.size();
    for (size_t i = 0; i < numRows; ++i)
    {
        outputSet.insert(theColumn[i]);
    }
}




/** \file ConvertTableColumnToSet.hpp
 * Converts a Table's column into a set.
 * \fn template<typename... T> void convertTableToSet(Table<T...>& theTable, set<typename tuple_element<0, Table<T...>::row_type>::type>& outputSet)
 * \section Usage Usage:
 * \code
 * IndexTable theTable;
 * // Add data to theTable;
 * set<size_t> theValues;
 * convertTableToSet(theTable, theValues);
 * \endcode
 */
template<typename... T>
void convertTableToSet(Table<T...>& theTable,
                       set<typename tuple_element<0, typename Table<T...>::row_type>::type>& outputSet)
{
#ifdef PROG_DEBUG
    if (tuple_size<typename Table<T...>::row_type>::value != 1)
    {
        cerr << "Warning.  convertTableToSet used on a table with multiple columns.  Consider using convertTableColumnToSet." << endl;
    }
#endif
    convertTableColumnToSet<0>(theTable, outputSet);
}



#endif

