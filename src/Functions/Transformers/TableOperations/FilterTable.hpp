#ifndef MOLBIOLIB_FILTERTABLE_H
#define MOLBIOLIB_FILTERTABLE_H 1

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
#include "src/Functions/Transformers/TableOperations/VectorDeleteTable.hpp"



/** \file FilterTable.hpp
 * Delete rows in a table.
 * \fn template<typename... T> void filterTable(Table<T...>& theTable, bool (*rowCriteria)(typename Table<T...>::row_type& theRow), bool fast = true)
 * Just delete rows based on a user-provided function.
 * \section Usage Usage:
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * myFastaTable myTable;
 * \endcode
 * Suppose one had the criteria function
 * \code
 * bool myCheckRows(myFastaTable::row_type& theRow1) { ... }
 * \endcode
 * then one can call #filterTable:
 * \code
 * filterTable(myTable, myCheckRows, false);
 * \endcode
 * where the false above is to use the slow (less memory) version of
 * #vectorDeleteTable.  By default, this parameter is set to true.
 */
template<typename... T>
void filterTable(Table<T...>& theTable,
                 bool (*rowCriteria)(typename Table<T...>::row_type& theRow),
                 bool fast = true)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_type;

    size_t numRows = theTable.size();
    vector<bool> delVector;
    delVector.reserve(numRows);

    row_type tempRow;
    for (size_t i = 0; i < numRows; ++i)
    {
        theTable.getRow(i, tempRow);
        delVector.push_back(rowCriteria(tempRow));
    }

    vectorDeleteTable(theTable, delVector, fast);
}




#endif

