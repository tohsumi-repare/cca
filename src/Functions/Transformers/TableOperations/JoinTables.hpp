#ifndef MOLBIOLIB_JOINTABLES_H
#define MOLBIOLIB_JOINTABLES_H 1

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
#include "src/Functions/Transformers/TupleOperations/JoinTuples.hpp"



template<typename T, typename U>
bool defaultInnerJoinTablesCompare(T& rowA, U& rowB)
{
    return (rowA == rowB);
}

template<typename T, typename U, typename V>
bool defaultInnerJoinTablesJoiner(T& rowA, U& rowB, V& rowOut)
{
    joinTuples(rowA, rowB, rowOut);
    return true;
}

template<typename T>
void defaultInnerJoinTablesEmptyFillIn(T& theRow)
{
    // Empty
}



/** \file JoinTables.hpp
 * Joins two #Table objects.
 * Inner joins two #Table objects.
 * \fn template<typename... T, typename... U, typename... V> void innerJoinTables(Table<T...>& tableA, Table<U...>& tableB, Table<V...>& tableOut, bool (*compare)(typename Table<T...>::row_type& rowA, typename Table<U...>::row_type& rowB) = defaultInnerJoinTablesCompare, bool (*joiner)(typename Table<T...>::row_type& rowA, typename Table<U...>::row_type& rowB, typename Table<V...>::row_type& rowC) = defaultInnerJoinTablesJoiner, void (*fillIn)(typename Table<V...>::row_type& theRow) = nullptr, size_t allocAmt = 0)
 * \section Usage Usage example:
 * \code
 * Table<string, string, int> table1, table2;
 * Table<string, string, int, string, string, int> table3;
 * \endcode
 * The more general version slower version of the form
 * \code
 * innerJoinTables<1,0>(table1, table2, table3, false, 9);
 * \endcode
 * is
 * \code
 * innerJoinTables(table1, table2, table3, comparator, joiner, fillIn, 9);
 * \endcode
 * One may provide the functions:
 * \code
 * bool comparator(Table<string, string, int>::row_type& rowA,
 *                 Table<string, string, int>::row_type& rowB) { ... }
 * bool joiner(Table<string, string, int>::row_type& rowA,
 *             Table<string, string, int>::row_type& rowB,
 *             Table<string, string, int, string, string, int>::row_type& rowOut) { ... }
 * \endcode
 * where the above is true if the resultant row is actually to be kept.
 * \code
 * void fillIn(Table<string, string, int, string, string, int>::row_type& theRow) { ... }
 * \endcode
 * where if not provided the default versions are called.  (See code.)  The
 * default behavior is where it must be table1 and table2 are of the same
 * type and a tuple comparison is done with a joinTuple done and unmatched
 * rows are not kept.  Unmatched rows are only kept if fillin is specified.
 * Use the following line, if you want the default but also want to keep the
 * unmatched lines.
 * \code
 * innerJoinTables(table1, table2, table3,
 *                 defaultInnerJoinTablesCompare,
 *                 defaultInnerJoinTablesJoiner,
 *                 defaultInnerJoinTablesEmptyFillIn);
 * \endcode
 */
template<typename... T, typename... U, typename... V>
void innerJoinTables(Table<T...>& tableA, Table<U...>& tableB,
                     Table<V...>& tableOut,
                     bool (*compare)(typename Table<T...>::row_type& rowA,
                                     typename Table<U...>::row_type& rowB) = defaultInnerJoinTablesCompare,
                     bool (*joiner)(typename Table<T...>::row_type& rowA,
                                    typename Table<U...>::row_type& rowB,
                                    typename Table<V...>::row_type& rowC) = defaultInnerJoinTablesJoiner,
                     void (*fillIn)(typename Table<V...>::row_type& theRow) = nullptr,
                     size_t allocAmt = 0)
{
#ifdef PROG_DEBUG
    assert(tableA.sameRowSizes() == true);
    assert(tableB.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_typeA;
    typedef typename Table<U...>::row_type row_typeB;
    typedef typename Table<V...>::row_type row_typeOut;

    size_t numRowsA = tableA.size();
    size_t numRowsB = tableB.size();

    tableOut.clear();
    if (allocAmt == 0)
    {
        tableOut.reserve(numRowsA * numRowsB);
    }
    else
    {
        tableOut.reserve(allocAmt);
    }


    for (size_t i = 0; i < tableA.labels.size(); ++i)
    {
        tableOut.labels[i] = tableA.labels[i];
    }
    for (size_t i = 0; i < tableB.labels.size(); ++i)
    {
        tableOut.labels[i+tableA.labels.size()] = tableB.labels[i];
    }


    row_typeA tempRowA;
    row_typeB tempRowB;
    row_typeOut tempRowOut;
    bool found = false, goodJoin = true;
    for (size_t i = 0; i < numRowsA; ++i)
    {
        tableA.getRow(i, tempRowA);
        found = false;
        for (size_t j = 0; j < numRowsB; ++j)
        {
            tableB.getRow(j, tempRowB);
            if (compare(tempRowA, tempRowB))
            {
                found = true;
                goodJoin = joiner(tempRowA, tempRowB, tempRowOut);
                if (goodJoin)
                {
                    tableOut.push_back(tempRowOut);
                }
            }
        }
        if (!found && fillIn != nullptr)
        {
            joinTuples(tempRowA, tempRowOut);
            fillIn(tempRowOut);
            tableOut.push_back(tempRowOut);
        }
    }

}





/** \file JoinTables.hpp
 * Joins two tables.
 * \fn template<size_t M, size_t N, typename... T, typename... U, typename... V> void innerJoinTables(Table<T...>& tableA, Table<U...>& tableB, Table<V...>& tableOut, bool keepNoMatches = false, size_t allocAmt = 0)
 * Inner joins two #Table objects.
 * \section Usage Usage example:
 * \code
 * Table<string, string, int> table1, table2;
 * Table<string, string, int, string, string, int> table3;
 * \endcode
 * The below does an inner join on key in the 2nd and 1st column of
 * table1 and table2 respectively and outputs the result to table3.
 * The false (default value) denotes that lines with no matches are deleted.
 * If true, then the lines are kept.  The below calls a faster (though using
 * more memory) version than the more general version.
 * \code
 * innerJoinTables<1,0>(table1, table2, table3, false, 9);
 * \endcode
 * where, 9, an option, is the number of entries to reserve in table 3.
 */
template<size_t M, size_t N, typename... T, typename... U, typename... V>
void innerJoinTables(Table<T...>& tableA, Table<U...>& tableB,
                     Table<V...>& tableOut, bool keepNoMatches = false,
                     size_t allocAmt = 0)
{
#ifdef PROG_DEBUG
    assert(tableA.sameRowSizes() == true);
    assert(tableB.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_typeA;
    typedef typename tuple_element<M,row_typeA>::type key_type;
    typedef typename Table<U...>::row_type row_typeB;
    typedef typename Table<V...>::row_type row_typeOut;

    size_t numRowsA = tableA.size();
    size_t numRowsB = tableB.size();

    tableOut.clear();
    if (allocAmt == 0)
    {
        tableOut.reserve(numRowsA * numRowsB);
    }
    else
    {
        tableOut.reserve(allocAmt);
    }

    for (size_t i = 0; i < tableA.labels.size(); ++i)
    {
        tableOut.labels[i] = tableA.labels[i];
    }
    for (size_t i = 0; i < tableB.labels.size(); ++i)
    {
        tableOut.labels[i+tableA.labels.size()] = tableB.labels[i];
    }

    row_typeA tempRowA;
    row_typeB tempRowB;
    row_typeOut tempRowOut;

    multimap<key_type, row_typeB> theMapB;
    for (size_t i = 0; i < numRowsB; ++i)
    {
        tableB.getRow(i, tempRowB);
        theMapB.insert(pair<key_type, row_typeB>(get<N>(tempRowB), tempRowB));
    }


    pair<typename multimap<key_type, row_typeB>::iterator,
         typename multimap<key_type, row_typeB>::iterator> theHits;
    typename multimap<key_type, row_typeB>::iterator iter;
    for (size_t i = 0; i < numRowsA; ++i)
    {
        tableA.getRow(i, tempRowA);
        if (theMapB.find(get<M>(tempRowA)) != theMapB.end())
        {
            theHits = theMapB.equal_range(get<M>(tempRowA));
            for (iter = theHits.first; iter != theHits.second; ++iter)
            {
                joinTuples(tempRowA, iter->second, tempRowOut);
                tableOut.push_back(tempRowOut);
            }
        }
        else if (keepNoMatches)
        {
            joinTuples(tempRowA, tempRowOut);
            tableOut.push_back(tempRowOut);
        }
    }
}




/** \file JoinTables.hpp
 * Joins two tables.
 * \fn template<typename... T, typename... U, typename... V> void outerJoinTables(Table<T...>& tableA, Table<U...>& tableB, Table<V...>& tableOut)
 * Outer join on two tables.
 * The below does an outer join on table1 and table2 and store in table3.
 * An outer join is where we make a copy of table2 for each row of table1
 * and then concatenate the columns into table3.
 * \section Usage Usage example:
 * With
 * \code
 * Table<string, string, int> table1, table2;
 * Table<string, string, int, string, string, int> table3;
 * \endcode
 * call
 * \code
 * outerJoinTables(table1, table2, table3);
 * \endcode
 */
template<typename... T, typename... U, typename... V>
void outerJoinTables(Table<T...>& tableA, Table<U...>& tableB,
                     Table<V...>& tableOut)
{
#ifdef PROG_DEBUG
    assert(tableA.sameRowSizes() == true);
    assert(tableB.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_typeA;
    typedef typename Table<U...>::row_type row_typeB;
    typedef typename Table<V...>::row_type row_typeOut;

    size_t numRowsA = tableA.size();
    size_t numRowsB = tableB.size();

    tableOut.clear();
    tableOut.reserve(numRowsA * numRowsB);

    for (size_t i = 0; i < tableA.labels.size(); ++i)
    {
        tableOut.labels[i] = tableA.labels[i];
    }
    for (size_t i = 0; i < tableB.labels.size(); ++i)
    {
        tableOut.labels[i+tableA.labels.size()] = tableB.labels[i];
    }

    row_typeA tempRowA;
    row_typeB tempRowB;
    row_typeOut tempRowOut;
    for (size_t i = 0; i < numRowsA; ++i)
    {
        tableA.getRow(i, tempRowA);
        for (size_t j = 0; j < numRowsB; ++j)
        {
            tableB.getRow(j, tempRowB);
            joinTuples(tempRowA, tempRowB, tempRowOut);
            tableOut.push_back(tempRowOut);
        }
    }
}


#endif

