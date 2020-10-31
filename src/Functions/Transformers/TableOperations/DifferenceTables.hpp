#ifndef MOLBIOLIB_DIFFERENCETABLES_H
#define MOLBIOLIB_DIFFERENCETABLES_H 1

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


/** \file DifferenceTables.hpp
 * Differences of two #Table objects.
 * \fn template<typename... T> void differenceTables(Table<T...>& tableA, Table<T...>& tableB, Table<T...>& tableC, bool delAllCopies = true, bool fast = true)
 * This returns a table that is the difference between two tables, i.e.
 * given tables A and B, both of the same type, compute A - B and return
 * in a new #Table (where the entries are sorted in some tuple ordering.
 * <b>Warning:</b> these operations, particularly
 * differenceTable, takes up a lot
 * of memory - a little more than the amount of memory to store the
 * actual tables.  This is also a slow process, O(n*log(n)), where
 * n is the number of rows.
 * Note that union is given by concatenateTableRows and
 * intersection is given by innerJoinTables.
 * \section Usage Usage:
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * myFastaTable table1, table2, tableOut;
 * ...
 * differenceTables(table1, table2, tableOut, true, false);
 *   // Compute table1 - table2 and store in tableOut.
 *   // The true above means that suppose table1 has 3 identical rows that
 *   // appear only once in table2, then all three rows are
 *   // deleted.  If false, then only one row in table1 is deleted.  (So,
 *   // three identical rows in table2 would be needed to
 *   // get rid of all same rows in table1.)  Default is true.
 *   // The false above means to use the slower version (but less memory) of
 *   // vectorDeleteTable.  By default, true (faster, but uses more memory).
 * \endcode
 * \section Implementation Implementation:
 * C = A - B
 */
template<typename... T>
void differenceTables(Table<T...>& tableA, Table<T...>& tableB,
                      Table<T...>& tableC,
                      bool delAllCopies = true, bool fast = true)
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
    vector<bool> delVector(numRows, false);
    for (size_t i = 0; i < numRows; ++i)
    {
        tableA.getRow(i, tempRow);
        if (tableBSet.find(tempRow) != tableBSet.end())
        {
            if (!delAllCopies)
            {
                iter = tableBSet.find(tempRow);
                tableBSet.erase(iter);
            }
            delVector[i] = true;
        }
    }

    tableC.clear();
    tableC = tableA;

    vectorDeleteTable(tableC, delVector, fast);
}





/** \file DifferenceTables.hpp
 * Differences of two #Table objects.
 * \fn template<size_t N, typename... T> void differenceTables(Table<T...>& tableA, Table<T...>& tableB, Table<T...>& tableC, bool delAllCopies = true, bool fast = true)
 * This returns a table that is the difference between two tables, i.e.
 * given tables A and B, both of the same type, compute A - B and return
 * in a new #Table (where the entries are sorted in some tuple ordering.
 * <b>Warning:</b> these
 * operations, particularly differenceTable, takes up a lot
 * of memory - a little more than the amount of memory to store the
 * actual tables.  This is also a slow process, O(n*log(n)), where
 * n is the number of rows.
 * Note that union is given by concatenateTableRows and
 * intersection is given by innerJoinTables.
 * \section Usage Usage:
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * myFastaTable table1, table2, tableOut;
 * ...
 * differenceTables<2>(table1, table2, tableOut);  // compute table1 - table2
 *                                                 // where rows are compared
 *                                                 // just in the third column.
 * \section Implementation Implementation:
 * C = A - B
 * \endcode
 */
template<size_t N, typename... T>
void differenceTables(Table<T...>& tableA, Table<T...>& tableB,
                      Table<T...>& tableC,
                      bool delAllCopies = true, bool fast = true)
{
#ifdef PROG_DEBUG
    assert(tableA.sameRowSizes() == true);
    assert(tableB.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_type;
    typedef typename tuple_element<N,row_type>::type key_type;
    multiset<key_type> tableBSet;
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
    vector<bool> delVector(numRows, false);
    for (size_t i = 0; i < numRows; ++i)
    {
        tableA.getRow(i, tempRow);
        if (tableBSet.find(get<N>(tempRow)) != tableBSet.end())
        {
            if (!delAllCopies)
            {
                iter = tableBSet.find(get<N>(tempRow));
                tableBSet.erase(iter);
            }
            delVector[i] = true;
        }
    }

    tableC.clear();
    tableC = tableA;

    vectorDeleteTable(tableC, delVector, fast);
}






/** \file DifferenceTables.hpp
 * Differences of two #Table objects.
 * \fn template<typename... T, typename... U> void differenceTables(Table<T...>& tableA, Table<U...>& tableB, Table<T...>& tableC, bool (*comparitor)(typename Table<T...>::row_type& rowA, typename Table<U...>::row_type& rowB), bool delAllCopies = true, bool fast = true)
 * This returns a table that is the difference between two tables, i.e.
 * given tables A and B, both of the same type, compute A - B and return
 * in a new #Table (where the entries are sorted in some tuple ordering.
 * <b>Warning:</b> these
 * operations, particularly differenceTable, takes up a lot
 * of memory - a little more than the amount of memory to store the
 * actual tables.  This is also a slow process, O(n*log(n)), where
 * n is the number of rows.
 * Note that union is given by concatenateTableRows and
 * intersection is given by innerJoinTables.
 * \section Usage Usage:
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * typedef Table<string, string, unsigned long, string> myNewTable;
 * myFastaTable table1, table2, tableOut;
 * \endcode
 * if one has a compare function:
 * \code
 * bool myCompare(myFastaTable::row_type& rowA,
 *                myNewTable::row_type& rowB) { .... }
 * \endcode
 * then one can do
 * \code
 * differenceTables(table1, table2, tableOut, myCompare);
 *    // where the above uses the comparator to compare tables of different
 *    // types and output the rows of table1 in tableOut such that the
 *    // comparator turns out false for all rows of table2.  Note that since
 *    // no ordering can be imposed between tables, this version of
 *    // differenceTables is very slow, O(n^2), and therefore should only be
 *    // used on smaller tables.  Also, this uses much more memory, as it needs
 *    // a copy of both tables kept.
 * \endcode
 */
template<typename... T, typename... U>
void differenceTables(Table<T...>& tableA, Table<U...>& tableB,
                      Table<T...>& tableC,
                      bool (*comparitor)(typename Table<T...>::row_type& rowA,
                              typename Table<U...>::row_type& rowB),
                      bool delAllCopies = true, bool fast = true)
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

    vector<bool> delVector(numRowsA, false);
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
                if (!delAllCopies)
                {
                    tempTableB.erase(j);
                }
                delVector[i] = true;
                found = true;
            }
        }
    }

    tableC.clear();
    tableC = tableA;

    vectorDeleteTable(tableC, delVector, fast);
}


#endif

