#ifndef MOLBIOLIB_VECTORDELETETABLE_H
#define MOLBIOLIB_VECTORDELETETABLE_H 1

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
#include "src/Functions/Algorithms/Tables/SortTable.hpp"



template<typename... T>
void vectorDeleteTableFast(IndexTable& delVector,
                           Table<T...>& theTable)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    if (delVector.size() == 0)
    {
        return;
    }
    typedef typename Table<T...>::row_type row_type;
    sortTable<0>(delVector);

    size_t numToDel = delVector.size();

    Table<T...> tempTable;

    for (size_t i = 0; i < tuple_size<row_type>::value; ++i)
    {
        tempTable.labels[i] = theTable.labels[i];
    }

    size_t tableSize = theTable.size();
    tempTable.resize(tableSize - numToDel);

    size_t currentDel = 0, currentTemp = 0;
    row_type tempRow;
    for (size_t i = 0; i < tableSize; ++i)
    {
        // First part of OR because means we have no more to delete...
        // Using short-circuit evaluation so should have no seg. fault.
        if (currentDel == numToDel || i != delVector[currentDel])
        {
            theTable.getRow(i, tempRow);
            tempTable.setRow(currentTemp, tempRow);
            ++currentTemp;
        }
        else
        {
            ++currentDel;
        }
    }

    theTable.clear();
    theTable = tempTable;
    tempTable.clear();
}



template<typename... T>
void vectorDeleteTableSlow(IndexTable& delVector, Table<T...>& theTable)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    if (delVector.size() == 0)
    {
        return;
    }
    sortTable<0>(delVector, true);

    // Want to delete backwards since only later entries affected.
    size_t numEntries = delVector.size();
    for (size_t i = 0; i < numEntries; ++i)
    {
        theTable.erase(delVector[i]);
    }
}



template<typename... T>
void vectorDeleteTable(Table<T...>& theTable, IndexTable& delVector,
                       bool fast = true)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    if (fast)
    {
        vectorDeleteTableFast(delVector, theTable);
    }
    else
    {
        vectorDeleteTableSlow(delVector, theTable);
    }
}




/** \file VectorDeleteTable.hpp
 * Delete many entries in a #Table.
 * \fn template<typename... T> void vectorDeleteTable(Table<T...>& theTable, vector<bool>& delVector, bool fast = true)
 * This routine does many deletions within a #Table given a vector of
 * indices to be deleted.  There are two modes.  One is fast and
 * does require memory to store a second table, and one is much slower
 * and requires no additional memory.
 * \section Usage Usage:
 * Given
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * myFastaTable myTable;
 * IndexTable entriesToBeDeleted;
 * // IndexTable is Table<size_t> typedef'd.
 * \endcode
 * where the above has the indices of which rows to be deleted
 * or alternatively
 * \code
 * Table<bool> entriesToBeDeleted;
 * \endcode
 * where the above is a vector of length number of rows
 * in #Table where if true, then delete that row.
 * Then, to delete entries in <code>myTable</code>, do
 * \code
 * vectorDeleteTable(myTable, entriesToBeDeleted, false);
 * \endcode
 * where the <code>false</code> specifies to use the
 * slow mode.  Default value is true.
 */
template<typename... T>
void vectorDeleteTable(Table<T...>& theTable, vector<bool>& delVector,
                       bool fast = true)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    size_t tableSize = theTable.size();
    Table<size_t> delSizeTypeVector;
    delSizeTypeVector.reserve(tableSize);
    for (size_t i = 0; i < tableSize; ++i)
        if (delVector[i])
        {
            delSizeTypeVector.push_back(i);
        }
    vectorDeleteTable(theTable, delSizeTypeVector, fast);
}






#endif

