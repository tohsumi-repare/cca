#ifndef MOLBIOLIB_CONCATENATETABLES_H
#define MOLBIOLIB_CONCATENATETABLES_H 1

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


// concatenateTablesRows(table2, table1) has been made obsolete by the
// Table's push_back(&Table theTable) function.



template<size_t N, size_t M, typename C, typename D>
class ConcatTablesCols
{
public:
    static void concat(C& vectuples1, D& vectuples2)
    {
        size_t numRows = get<N>(vectuples1).size();
        get<N+M>(vectuples2).clear();
        get<N+M>(vectuples2).reserve(numRows);
        for (size_t i = 0; i < numRows; ++i)
        {
            get<N+M>(vectuples2).push_back(get<N>(vectuples1)[i]);
        }
        ConcatTablesCols<N-1,M,C,D>::concat(vectuples1, vectuples2);
    }
};
template<size_t M, typename C, typename D>
class ConcatTablesCols<0,M,C,D>
{
public:
    static void concat(C& vectuples1, D& vectuples2)
    {
        size_t numRows = get<0>(vectuples1).size();
        get<M>(vectuples2).clear();
        get<M>(vectuples2).reserve(numRows);
        for (size_t i = 0; i < numRows; ++i)
        {
            get<M>(vectuples2).push_back(get<0>(vectuples1)[i]);
        }
    }
};



/** \file ConcatenateTables.hpp
 * Add rows to end of one Table from another Table.
 * \fn template<typename... T, typename... U, typename... V> void concatenateTablesCols(Table<T...>& table1, Table<U...>& table2, Table<V...>& table3)
 * \section Usage Usage example:
 * \code
 * Table<string, string, int> table1, table2, table5;
 * Table<string, double> table3;
 * Table<string, double, string, string, int> table4;
 * \endcode
 * Tacks together the columns of table3 and table1 outputs it into table4.
 * \code
 * concatenateTablesCols(table3, table1, table4);
 * \endcode
 * \section ProgNotes Programming Notes:
 * Did not use joinTuples, since the below implementation
 * was written before joinTuples was and is faster anyway.
 */
template<typename... T, typename... U, typename... V>
void concatenateTablesCols(Table<T...>& table1, Table<U...>& table2,
                           Table<V...>& table3)
{
    typedef typename Table<T...>::row_type tableRow_typeFrom1;
    typedef typename Table<T...>::cols_type tableCols_typeFrom1;
    typedef typename Table<U...>::row_type tableRow_typeFrom2;
    typedef typename Table<U...>::cols_type tableCols_typeFrom2;
    typedef typename Table<V...>::cols_type tableCols_typeTo;

    for (size_t i = 0; i < tuple_size<tableRow_typeFrom1>::value; ++i)
    {
        table3.labels[i] = table1.labels[i];
    }
    for (size_t i = 0; i < tuple_size<tableRow_typeFrom2>::value; ++i)
    {
        table3.labels[i+tuple_size<tableRow_typeFrom1>::value] = table2.labels[i];
    }

    ConcatTablesCols<tuple_size<tableRow_typeFrom1>::value - 1,
                     0,
                     tableCols_typeFrom1,
                     tableCols_typeTo>::concat(table1.theData,
                             table3.theData);

    ConcatTablesCols<tuple_size<tableRow_typeFrom2>::value - 1,
                     tuple_size<tableRow_typeFrom1>::value,
                     tableCols_typeFrom2,
                     tableCols_typeTo>::concat(table2.theData,
                             table3.theData);
}


#endif

