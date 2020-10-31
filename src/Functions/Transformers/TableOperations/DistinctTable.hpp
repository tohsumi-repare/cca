#ifndef MOLBIOLIB_DISTINCTTABLE_H
#define MOLBIOLIB_DISTINCTTABLE_H 1

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


// Not an unordered map since we do not know what type of key is being hashed.




/** \file DistinctTable.hpp
 * Make rows in #Table distinct #Table
 * \fn template<size_t N, typename... T> void distinctTable(Table<T...>& theTable, Table<typename tuple_element<N, typename Table<T...>::row_type>::type, size_t>& theFrequency, void (*rowMerger)(typename Table<T...>::row_type& r_source, typename Table<T...>::row_type& r_target) = nullptr)
 * A key from a column is taken from the table and rows are make distinct on
 * that.  Since the merge of two rows are required (perhaps one will wish to
 * use one of the columns as a count) but is an undefined process, a
 * user-written function is required.  After rows are merged, all the
 * non-distinct rows are deleted.  Thus, the table afterwards will typically
 * have less rows than the table before this is called.  The table is also
 * sorted by the key.  Also, outputs a #Table of key with frequency.  This
 * can be merged with the output using #innerJoinTables
 * \section Usage Usage:
 * Starting with
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * myFastaTable myTable;
 * \endcode
 * and the optional user-written merging function
 * \code
 * void myMergeRow(myFastaTable::row_type& r_source,
 *                 myFastaTable::row_type& r_target) {
 *    // merge contents of r_source into r_target
 * }
 * \endcode
 * then make <code>myTable</code> distinct using
 * \code
 * Table<string, size_t> myTableFrequency;
 * distinctTable<1>(myTable, myTableFrequency, myMergeRow);
 *                                           // Make distinct on keys in
 *                                           // the 2nd column.
 *    // Only works if all the columns have the same number of rows.
 *
 * Optionally, one may omit the function <code>myMergeRow</code> so that
 * during a merge, the new row to be merged will simply not be
 * used (retaining only the first row encountered).
 * \endcode
 */
template<size_t N, typename... T>
void distinctTable(Table<T...>& theTable,
                   Table<typename tuple_element<N, typename Table<T...>::row_type>::type,
                   size_t>& theFrequency,
                   void (*rowMerger)(typename Table<T...>::row_type& r_source,
                                     typename Table<T...>::row_type& r_target) = nullptr)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_type;
    typedef typename tuple_element<N, row_type>::type key_type;

    map<key_type, row_type> theMap;
    typedef typename map<key_type, row_type>::iterator map_iter;
    map<key_type, size_t> theFreq;
    typedef typename map<key_type, size_t>::iterator key_iter;


    // Go through rows, and populate theMap.  Then clear the table and
    // repopulate the table by iterating through theMap.

    size_t numRows = get<N>(theTable.theData).size();
    for (size_t i = 0; i < numRows; ++i)
    {
        row_type tempRow;
        theTable.getRow(i, tempRow);
        if (theMap.find(get<N>(tempRow)) == theMap.end())
        {
            theMap[get<N>(tempRow)] = tempRow;
            theFreq[get<N>(tempRow)] = 1;
        }
        else
        {
            if (rowMerger != nullptr)
            {
                rowMerger(tempRow, theMap[get<N>(tempRow)]);
            }
            theFreq[get<N>(tempRow)] += 1;
        }
    }

    theTable.clear();
    theTable.resize(theMap.size());

    size_t currentRow = 0;
    for (map_iter iter = theMap.begin(); iter != theMap.end(); ++iter)
    {
        theTable.setRow(currentRow, iter->second);
        ++currentRow;
    }


    theFrequency.clear();
    theFrequency.resize(theMap.size());

    currentRow = 0;
    for (key_iter iter = theFreq.begin(); iter != theFreq.end(); ++iter)
    {
        tuple<key_type, size_t> tempRow(iter->first, iter->second);
        theFrequency.setRow(currentRow, tempRow);
        ++currentRow;
    }
}





/** \file DistinctTable.hpp
 * Make rows in #Table distinct
 * \fn template<size_t N, typename... T> void distinctTable(Table<T...>& theTable, IndexTable& theFrequency, void (*rowMerger)(typename Table<T...>::row_type& r_source, typename Table<T...>::row_type& r_target) = nullptr)
 * Same as above except the frequency table is just the frequencies in the same
 * order as the table entries, but no key associated with it.  Thus, just pass
 * a #IndexTable for the frequency.
 */
template<size_t N, typename... T>
void distinctTable(Table<T...>& theTable,
                   IndexTable& theFrequency,
                   void (*rowMerger)(typename Table<T...>::row_type& r_source,
                                     typename Table<T...>::row_type& r_target) = nullptr)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_type;
    typedef typename tuple_element<N, row_type>::type key_type;

    map<key_type, row_type> theMap;
    typedef typename map<key_type, row_type>::iterator map_iter;
    map<key_type, size_t> theFreq;
    typedef typename map<key_type, size_t>::iterator key_iter;


    // Go through rows, and populate theMap.  Then clear the table and
    // repopulate the table by iterating through theMap.

    size_t numRows = get<N>(theTable.theData).size();
    for (size_t i = 0; i < numRows; ++i)
    {
        row_type tempRow;
        theTable.getRow(i, tempRow);
        if (theMap.find(get<N>(tempRow)) == theMap.end())
        {
            theMap[get<N>(tempRow)] = tempRow;
            theFreq[get<N>(tempRow)] = 1;
        }
        else
        {
            if (rowMerger != nullptr)
            {
                rowMerger(tempRow, theMap[get<N>(tempRow)]);
            }
            theFreq[get<N>(tempRow)] += 1;
        }
    }

    theTable.clear();
    theTable.resize(theMap.size());

    size_t currentRow = 0;
    for (map_iter iter = theMap.begin(); iter != theMap.end(); ++iter)
    {
        theTable.setRow(currentRow, iter->second);
        ++currentRow;
    }


    theFrequency.clear();
    theFrequency.resize(theMap.size());

    currentRow = 0;
    for (key_iter iter = theFreq.begin(); iter != theFreq.end(); ++iter)
    {
        tuple<size_t> tempRow(iter->second);
        theFrequency.setRow(currentRow, tempRow);
        ++currentRow;
    }
}





/** \file DistinctTable.hpp
 * Make rows in #Table distinct
 * \fn template<size_t N, typename... T> void distinctTable(Table<T...>& theTable, void (*rowMerger)(typename Table<T...>::row_type& r_source, typename Table<T...>::row_type& r_target) = nullptr)
 * This is the same as the above function, except without retaining the
 * frequency information.  The merge function is optional here as well.
 */
template<size_t N, typename... T>
void distinctTable(Table<T...>& theTable,
                   void (*rowMerger)(typename Table<T...>::row_type& r_source,
                                     typename Table<T...>::row_type& r_target) = nullptr)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_type;
    typedef typename tuple_element<N, row_type>::type key_type;

    map<key_type, row_type> theMap;
    typedef typename map<key_type, row_type>::iterator map_iter;


    // Go through rows, and populate theMap.  Then clear the table and
    // repopulate the table by iterating through theMap.

    size_t numRows = get<N>(theTable.theData).size();
    for (size_t i = 0; i < numRows; ++i)
    {
        row_type tempRow;
        theTable.getRow(i, tempRow);
        if (theMap.find(get<N>(tempRow)) == theMap.end())
        {
            theMap[get<N>(tempRow)] = tempRow;
        }
        else if (rowMerger != nullptr)
        {
            rowMerger(tempRow, theMap[get<N>(tempRow)]);
        }
    }

    theTable.clear();
    theTable.resize(theMap.size());

    size_t currentRow = 0;
    for (map_iter iter = theMap.begin(); iter != theMap.end(); ++iter)
    {
        theTable.setRow(currentRow, iter->second);
        ++currentRow;
    }

}


#endif

