#ifndef MOLBIOLIB_VECTORBOUNDS_H
#define MOLBIOLIB_VECTORBOUNDS_H 1

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


#include <vector>


/** \file VectorBounds.hpp
 * MolBioLib C++ analogs of lower_bound and upper_bound for vector.
 * \fn template<typename T> size_t lower_boundIndex(vector<T>& theVec, T value, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool (*compare)(T& x, T& y) = nullptr)
 * Finds the lower bound index on a vector.  This is the index for the first
 * element in the sorted range [first,last) which does not compare less than
 * value.  This takes the N'th column (0-index) of the table as the vector.
 * \section Usage Usage:
 * \code
 * vector<long> myVec;
 * ...
 * size_t low = lower_boundIndex<2>(myVec, 55);
 * \endcode
 * If the vector size is zero, numeric_limits<size_t>::max() is returned.  One
 * must check that the first element is indeed >= value, as 0 is the lowest
 * return value, if 0 is returned.
 * Optionally, one may pass in a comparison function which compares the values.
 * The comparison is equivalent to operator<.
 */
template<typename T>
size_t lower_boundIndex(vector<T>& theVec,
                        T value,
                        size_t first = 0,
                        size_t last = numeric_limits<size_t>::max(),
                        bool (*compare)(T& x, T& y) = nullptr)
{
    size_t result = numeric_limits<size_t>::max();
    size_t vecLength = theVec.size();
    if (vecLength == 0)
    {
        return result;    // No entries in table.
    }
    else
    {
        result = first;
    }
    if (last == numeric_limits<size_t>::max())
    {
        last = theVec.size();
    }
    size_t count = last - first;
    size_t check = 0, step = 0;
    while (count > 0)
    {
        check = result;
        step = count/2;
        check += step;
        if ((check < vecLength) &&
                (((compare == nullptr) && (theVec[check] < value)) ||
                 ((compare != nullptr) && (compare(theVec[check], value)))))
        {
            ++check;
            result = check;
            if (count > (step + 1))
            {
                count -= step + 1;
            }
            else
            {
                count = 0;
            }
        }
        else
        {
            count = step;
        }
    }
    return result;
}



/** \file VectorBounds.hpp
 * MolBioLib C++ analogs of lower_bound and upper_bound for vector.
 * \fn template<typename T> size_t upper_boundIndex(vector<T>& theVec, T value, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool (*compare)(T& x, T& y) = nullptr)
 * Finds the upper bound index on a vector.  This is the index for the first
 * element in the sorted range [first,last) which compares greater than value.
 * One must check that the first element is indeed > value, as 0 is the lowest
 * return value, if 0 is returned.
 * This takes the N'th column (0-index) of the table as the vector.
 * \section Usage Usage:
 * \code
 * vector<long> myVec;
 * ...
 * size_t up = upper_boundIndex<2>(myVec, 55);
 * \endcode
 * If the vector size is zero, numeric_limits<size_t>::max() is returned.
 * Optionally, one may pass in a comparison function which compares the values.
 * The comparison is equivalent to operator<.
 */
template<typename T>
size_t upper_boundIndex(vector<T>& theVec,
                        T value,
                        size_t first = 0,
                        size_t last = numeric_limits<size_t>::max(),
                        bool (*compare)(T& x, T& y) = nullptr)
{
    size_t result = numeric_limits<size_t>::max();
    size_t vecLength = theVec.size();
    if (vecLength == 0)
    {
        return result;    // No entries in table.
    }
    else
    {
        result = first;
    }
    if (last == numeric_limits<size_t>::max())
    {
        last = theVec.size();
    }
    size_t count = last - first;
    size_t check = 0, step = 0;
    while (count > 0)
    {
        check = result;
        step = count/2;
        check += step;
        if ((check < vecLength) &&
                (((compare == nullptr) && (!(value < theVec[check]))) ||
                 ((compare != nullptr) && (!(compare(value, theVec[check]))))))
        {
            ++check;
            result = check;
            if (count > (step + 1))
            {
                count -= step + 1;
            }
            else
            {
                count = 0;
            }
        }
        else
        {
            count = step;
        }
    }
    return result;
}


#endif

