#ifndef MOLBIOLIB_CHECKCONSTANTTABLE_H
#define MOLBIOLIB_CHECKCONSTANTTABLE_H 1

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


/** \file CheckConstantTable.hpp
 * Checks if a Table's column has the same value
 * by comparing the first element to all other elements in the column.
 * \fn template<size_t N, typename... T> bool checkConstantTable(Table<T...>& theTable, size_t begin = 0, size_t end = numeric_limits<size_t>::max())
 * Check if table column (casted to double) is numerically constant.
 * Return true if constant.
 */
template<size_t N, typename... T>
bool checkConstantTable(Table<T...>& theTable,
                        size_t begin = 0,
                        size_t end = numeric_limits<size_t>::max())
{
    if (end == numeric_limits<size_t>::max())
    {
        end = theTable.size() - 1;
    }

    if (begin == end)
    {
        return true;
    }
#ifdef DEBUG
    assert(begin < theTable.size());
    assert(begin <= end);
    assert(end < theTable.size());
#endif
    double first = static_cast<double>(get<N>(theTable.theData)[begin]);
    double temp = 0.0;
    for (size_t i = begin + 1; i <= end; ++i)
    {
        temp = static_cast<double>(get<N>(theTable.theData)[i]);
        if (temp != first)
        {
            return false;
        }
    }
    return true;

}


#endif

