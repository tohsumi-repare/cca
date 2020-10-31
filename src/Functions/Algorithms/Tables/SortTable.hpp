#ifndef MOLBIOLIB_SORTTABLE_H
#define MOLBIOLIB_SORTTABLE_H 1

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





template<size_t N, typename C>
class TableSwapRows
{
public:
    static void swapRows(C& theData, size_t& i, size_t& j)
    {
        TableSwapRows<N-1,C>::swapRows(theData, i, j);
        swap(get<N>(theData)[i], get<N>(theData)[j]);
    }
};
template<typename C>
class TableSwapRows<0,C>
{
public:
    static void swapRows(C& theData, size_t& i, size_t& j)
    {
        swap(get<0>(theData)[i], get<0>(theData)[j]);
    }
};



template<size_t N, typename... T>
class TableSingleColCompareRows
{
public:
    Table<T...>& table;
    TableSingleColCompareRows(Table<T...>& theTable) : table(theTable) {}
    bool operator() (size_t i, size_t j)
    {
        typedef typename Table<T...>::row_type row_type;
        row_type rowI, rowJ;
        table.getRow(i, rowI);
        table.getRow(j, rowJ);
        return (get<N>(rowI) < get<N>(rowJ));
    }
};




/** \file SortTable.hpp
 * Sort #Table.
 * \fn template<size_t N, typename... T> void sortTable(Table<T...>& theTable, bool reverse = false, bool stable = false)
 * Sorts the input table based on a comparator of two rows.
 * \section Usage Usage:
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * myFastaTable myTable;
 * ...
 * sortTable<2>(myTable, true, true);  // This will sort the table using the
 *                                     // third key (unsigned long) and will
 *                                     // reverse it (so the table will be in
 *                                     // descending order) afterwards.  By
 *                                     // default, reverse is false.  The last
 *                                     // true will do a stable sort.
 * \endcode
 */
template<size_t N, typename... T>
void sortTable(Table<T...>& theTable, bool reverse = false, bool stable = false)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_type;
    typedef typename Table<T...>::cols_type cols_type;

    size_t numRows = theTable.size();
    vector<size_t> tempPermutation(numRows);
    for (size_t i = 0; i < numRows; ++i)
    {
        tempPermutation[i] = i;
    }

    if (!stable)
    {
        sort(tempPermutation.begin(), tempPermutation.end(),
             TableSingleColCompareRows<N,T...>(theTable) );
    }
    else
    {
        stable_sort(tempPermutation.begin(), tempPermutation.end(),
                    TableSingleColCompareRows<N,T...>(theTable) );
    }

    vector<size_t> permutation(numRows);
    for (size_t i = 0; i < numRows; ++i)
    {
        permutation[tempPermutation[i]] = i;
    }

    tempPermutation.clear();


    for (size_t i = 0; i < numRows; ++i)
    {
        while (permutation[i] != i)
        {
            TableSwapRows<tuple_size<row_type>::value - 1,
                          cols_type>::swapRows(theTable.theData,
                                               i, permutation[i]);
            swap(permutation[i], permutation[permutation[i]]);
        }
    }

    if (reverse)
    {
        size_t targetRows = 0;
        if (numRows % 2 == 0)
        {
            targetRows = numRows/2;
        }
        else
        {
            targetRows = (numRows-1)/2;
        }
        for (size_t i = 0; i < targetRows; ++i)
        {
            // Forced to do a static cast below since while we know this is always
            // positive, the compiler does not, so the - turns this into an int.
            size_t j = static_cast<size_t>(numRows - 1 - i);
            TableSwapRows<tuple_size<row_type>::value - 1,
                          cols_type>::swapRows(theTable.theData, i, j);
        }
    }

    // permutation's deletion operator will be called at the end.

}








template<typename... T>
class TableCompareRows : public binary_function<size_t, size_t, bool>
{
public:
    Table<T...>& table;
    bool (*cmp)(typename Table<T...>::row_type& ri,
                typename Table<T...>::row_type& rj);
    TableCompareRows(Table<T...>& theTable,
                     bool (*comparator)(typename Table<T...>::row_type& ri,
                                        typename Table<T...>::row_type& rj)) :
        table(theTable), cmp(comparator) {}
    TableCompareRows() {}
    bool operator() (size_t i, size_t j)
    {
        typedef typename Table<T...>::row_type row_type;
        row_type rowI, rowJ;
        table.getRow(i, rowI);
        table.getRow(j, rowJ);
        return cmp(rowI, rowJ);
    }
};







/** \file SortTable.hpp
 * Sort #Table.
 * \fn template<typename... T> void sortTable(Table<T...>& theTable, bool (*cmp)(typename Table<T...>::row_type& rowA, typename Table<T...>::row_type& rowB), bool reverse = false, bool stable = false)
 * Sorts the input table based on a comparator of two rows.
 * \section Usage Usage:
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * myFastaTable myTable;
 * ...
 * \endcode
 * Suppose one had the comparative function
 * \code
 * bool myCompareRows(myFastaTable::row_type& r1,
 *                    myFastaTable::row_type& r2) { ... }
 * \endcode
 * where we know <code>row_type</code> is some tuple, such as for the
 * above example
 * <code>myFastaTable::row_type</code> can be replaced with
 * <code>tuple<string, string, unsigned long></code>, then we
 * could do a more general comparison by doing
 * \code
 * sortTable(myTable, myCompareRows);
 * \endcode
 * In the above case, if <code>myCompareRows</code> acts
 * like the <code><</code> operator, then the resultant table
 * will be in ascending order.
 */
template<typename... T>
void sortTable(Table<T...>& theTable,
               bool (*cmp)(typename Table<T...>::row_type& rowA,
                           typename Table<T...>::row_type& rowB),
               bool reverse = false, bool stable = false)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_type;
    typedef typename Table<T...>::cols_type cols_type;

    size_t numRows = theTable.size();
    vector<size_t> tempPermutation(numRows);
    for (size_t i = 0; i < numRows; ++i)
    {
        tempPermutation[i] = i;
    }
    if (!stable)
    {
        sort(tempPermutation.begin(), tempPermutation.end(),
             TableCompareRows<T...>(theTable, cmp) );
    }
    else
    {
        stable_sort(tempPermutation.begin(), tempPermutation.end(),
                    TableCompareRows<T...>(theTable, cmp) );
    }

    vector<size_t> permutation(numRows);
    for (size_t i = 0; i < numRows; ++i)
    {
        permutation[tempPermutation[i]] = i;
    }

    tempPermutation.clear();

    for (size_t i = 0; i < numRows; ++i)
    {
        while (permutation[i] != i)
        {
            TableSwapRows<tuple_size<row_type>::value - 1,
                          cols_type>::swapRows(theTable.theData,
                                               i, permutation[i]);
            swap(permutation[i], permutation[permutation[i]]);
        }
    }

    if (reverse)
    {
        size_t targetRows = 0;
        if (numRows % 2 == 0)
        {
            targetRows = numRows/2;
        }
        else
        {
            targetRows = (numRows-1)/2;
        }
        for (size_t i = 0; i < targetRows; ++i)
        {
            // Forced to do a static cast below since while we know this is always
            // positive, the compiler does not, so the - turns this into an int.
            size_t j = static_cast<size_t>(numRows - 1 - i);
            TableSwapRows<tuple_size<row_type>::value - 1,
                          cols_type>::swapRows(theTable.theData, i, j);
        }
    }

    // permutation's deletion operator will be called at the end.
}

#endif

