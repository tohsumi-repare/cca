#ifndef MOLBIOLIB_INTERSECTIONTABLES_H
#define MOLBIOLIB_INTERSECTIONTABLES_H 1

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



/** \file IntersectionTables.hpp
 * Intersection of two #Table objects.
 * \fn template<typename... T> void intersectionTables(Table<T...>& tableA, Table<T...>& tableB, Table<T...>& tableC, bool syncCopies = false, bool fast = true)
 * This returns a table that is the intersection between two tables, i.e.
 * givens tables A and B, both of the same type, compute A cap B and return
 * in a new table (where the entries are sorted in some tuple ordering.
 * <i>Warning: these operations takes up a lot of memory - a little more than
 * the amount of memory to store the actual tables.  This is also a
 * slow process, O(n*log(n)), where n is the number of rows.</i>
 * \section Usage Usage:
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * myFastaTable table1, table2, tableOut;
 * ...
 * intersectionTables(table1, table2, tableOut, false, false);
 *    // Compute table1 - table2 and store in tableOut.
 *    // The first false (default value) above means that if there is three
 *    // rows the same in table1 and just one row in table2 matching, then all
 *    // three rows are kept.  If the first false is instead replaced by true,
 *    // then the number of rows that gets put in the intersection
 *    // must match the lesser number of identical rows in table1 and table2.
 *    // The false above means to use the slower version (but less memory) of
 *    // vectorDeleteTable.  By default, true (faster, but uses more memory).
 * \endcode
 */
template<typename... T>
void intersectionTables(Table<T...>& tableA, Table<T...>& tableB,
                        Table<T...>& tableC,
                        bool syncCopies = false, bool fast = true)
{
#ifdef PROG_DEBUG
    assert(tableA.sameRowSizes() == true);
    assert(tableB.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_type;
    multiset<row_type> tableBSet;
    typename multiset<row_type>::iterator iter;
    size_t numRows = tableB.size();

    for (size_t i = 0; i < tableA.labels.size(); ++i)
    {
        tableC.labels[i] = tableA.labels[i];
    }

    row_type tempRow;
    for (size_t i = 0; i < numRows; ++i)
    {
        tableB.getRow(i, tempRow);
        tableBSet.insert(tempRow);
    }

    numRows = tableA.size();
    vector<bool> delVector(numRows, true);
    for (size_t i = 0; i < numRows; ++i)
    {
        tableA.getRow(i, tempRow);
        if (tableBSet.find(tempRow) != tableBSet.end())
        {
            if (syncCopies)
            {
                iter = tableBSet.find(tempRow);
                tableBSet.erase(iter);
            }
            delVector[i] = false;
        }
    }

    tableC.clear();
    tableC = tableA;

    vectorDeleteTable(tableC, delVector, fast);
}




/** \file IntersectionTables.hpp
 * Intersection of two #Table objects.
 * \fn template<size_t N, typename... T> void intersectionTables(Table<T...>& tableA, Table<T...>& tableB, Table<T...>& tableC, bool syncCopies = false, bool fast = true)
 * This returns a table that is the intersection between two tables, i.e.
 * givens tables A and B, both of the same type, compute A cap B and return
 * in a new table (where the entries are sorted in some tuple ordering.
 * <i>Warning: these operations takes up a lot of memory - a little more than
 * the amount of memory to store the actual tables.  This is also a
 * slow process, O(n*log(n)), where n is the number of rows.</i>
 * \section Usage Usage:
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * myFastaTable table1, table2, tableOut;
 * ...
 * intersectionTables<2>(table1, table2, tableOut);  // compute table1 - table2
 *                                                   // where rows are compared
 *                                                   // just in the third
 *                                                   // column.
 * \endcode
 */
template<size_t N, typename... T>
void intersectionTables(Table<T...>& tableA, Table<T...>& tableB,
                        Table<T...>& tableC,
                        bool syncCopies = false, bool fast = true)
{
#ifdef PROG_DEBUG
    assert(tableA.sameRowSizes() == true);
    assert(tableB.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_type;
    typedef typename tuple_element<N,row_type>::type key_type;
    set<key_type> tableBSet;
    typename multiset<key_type>::iterator iter;
    size_t numRows = tableB.size();

    for (size_t i = 0; i < tableA.labels.size(); ++i)
    {
        tableC.labels[i] = tableA.labels[i];
    }

    row_type tempRow;
    for (size_t i = 0; i < numRows; ++i)
    {
        tableB.getRow(i, tempRow);
        tableBSet.insert(get<N>(tempRow));
    }

    numRows = tableA.size();
    vector<bool> delVector(numRows, true);
    for (size_t i = 0; i < numRows; ++i)
    {
        tableA.getRow(i, tempRow);
        if (tableBSet.find(get<N>(tempRow)) != tableBSet.end())
        {
            if (syncCopies)
            {
                iter = tableBSet.find(get<N>(tempRow));
                tableBSet.erase(iter);
            }
            delVector[i] = false;
        }
    }

    tableC.clear();
    tableC = tableA;

    vectorDeleteTable(tableC, delVector, fast);
}






/** \file IntersectionTables.hpp
 * Intersection of two #Table objects.
 * \fn template<typename... T, typename... U> void intersectionTables(Table<T...>& tableA, Table<U...>& tableB, Table<T...>& tableC, bool (*comparitor)(typename Table<T...>::row_type& rowA, typename Table<U...>::row_type& rowB), bool syncCopies = false, bool fast = true)
 * This returns a table that is the intersection between two tables, i.e.
 * givens tables A and B, both of the same type, compute A cap B and return
 * in a new table (where the entries are sorted in some tuple ordering.
 * <i>Warning: these operations takes up a lot of memory - a little more than
 * the amount of memory to store the actual tables.  This is also a
 * slow process, O(n*log(n)), where n is the number of rows.</i>
 * \section Usage Usage:
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * myFastaTable table1, table2, tableOut;
 * ...
 * \endcode
 * if one has a compare function:
 * \code
 * bool myCompare(myFastaTable::row_type& rowA,
 *                myNewTable::row_type& rowB) { .... }
 * \endcode
 * then one can do
 * \code
 * intersectionTables(table1, table2, tableOut, myCompare);
 *    // where the above uses the comparator to compare tables of different
 *    // types and output the rows of table1 in tableOut such that the
 *    // comparator turns out true for all rows of table2.  Note that since
 *    // no ordering can be imposed between tables, this version of
 *    // intersectionTables is very slow, O(n^2), and therefore should only be
 *    // used on smaller tables.
 * \endcode
 */
template<typename... T, typename... U>
void intersectionTables(Table<T...>& tableA, Table<U...>& tableB,
                        Table<T...>& tableC,
                        bool (*comparitor)(typename Table<T...>::row_type& rowA,
                                typename Table<U...>::row_type& rowB),
                        bool syncCopies = false, bool fast = true)
{
#ifdef PROG_DEBUG
    assert(tableA.sameRowSizes() == true);
    assert(tableB.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_typeA;
    typedef typename Table<U...>::row_type row_typeB;

    for (size_t i = 0; i < tableA.labels.size(); ++i)
    {
        tableC.labels[i] = tableA.labels[i];
    }

    size_t numRowsA = tableA.size();
    size_t numRowsB = tableB.size();

    row_typeA tempRowA;
    row_typeB tempRowB;

    Table<U...> tempTableB = tableB;

    vector<bool> delVector(numRowsA, true);
    bool found = false;
    for (size_t i = 0; i < numRowsA; ++i)
    {
        tableA.getRow(i, tempRowA);
        found = false;
        for (size_t j = 0; j < numRowsB && !found; ++j)
        {
            tempTableB.getRow(j, tempRowB);
            if (comparitor(tempRowA, tempRowB))
            {
                if (syncCopies)
                {
                    tempTableB.erase(j);
                }
                delVector[i] = false;
                found = true;
            }
        }
    }

    tableC.clear();
    tableC = tableA;

    vectorDeleteTable(tableC, delVector, fast);
}


#endif

