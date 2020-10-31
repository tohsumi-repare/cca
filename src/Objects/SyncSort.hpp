#ifndef MOLBIOLIB_SYNCSORT_H
#define MOLBIOLIB_SYNCSORT_H 1

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


#include "src/include/PrimitiveTypes.hpp"


/** \file SyncSort.hpp
 * Sort a number of vectors based on the sorted ordering of one vector
 */


/** Contains methods and appliers
 * The below is one templated class.
 * Done as a class so as to hold the
 * permutation vector.
 *
 * \section Usage Usage:
 * \code
 * vector<int> keys;  // fill in...
 * vector<string> data1;  // fill in...
 * vector<double> data2;  // fill in...
 * SyncSort<int> keysSort(keys, true, false, true, true);
 * // true, false, true, and true are the default parameters, described below.
 * // The first "true" means to sort from small to largest.  If false, then
 * // sort in reverse order.
 * // The second "false" means not to distinctify on keys.  If true, then
 * // distinctify keys and delete same elements when applied to other vectors.
 * // The third "true" means to keep the first one.  If false, keep last one.
 * // No need to do below to keys - automatically done.
 * // The last "true" means to use a stable sort.  If false, use a non-stable
 * // sort which is faster.
 * keysSort.apply(data1);
 * keysSort.apply(data2);
 * keysSort.clear();   // If the scope of keysSort is not local
 *                     // and not using anymore.
 *
 * // Can also set the above parameters later with
 * // setForward(bool), setDistinct(bool), setKeepFirst(bool),
 * // and setStable(bool).
 *
 * // Can change keys and resync by doing
 * keysSort.sortKeys();
 * keysSort.apply(data1);
 * keysSort.apply(data2);
 * \endcode
 *
 * Note that keeping the first key in forward order is the same as keeping
 * the last key in reverse order and vice-versa.
 *
 * Programmer's note: when in case of having things like above where all
 * vectors are the same size, it is preferable to do
 * \code
 * Table<int, string, double> keysAndData;  // fill in...
 * sortTable<0>(keysAndData);
 * \endcode
 *
 */
template<typename T>
class SyncSort
{
public:
    template<typename V>
    void apply(vector<V>& theVec)
    {
#ifdef PROG_DEBUG
        assert(theVec.size() == numRows);
#endif
        // Below: do not permute if there is nothing to do
        //    Also, if numRows = 0, numRows-1 is numeric_limits<size_t>::max().
        if (numRows <= 1)
        {
            return;
        }
        vector<size_t> tempPermutation = permutationVector;
        for (size_t i = 0; i < numRows; ++i)
        {
            while (tempPermutation[i] != i)
            {
                swap(theVec[i], theVec[tempPermutation[i]]);
                swap(tempPermutation[i], tempPermutation[tempPermutation[i]]);
            }
        }

        if (makeDistinct)
        {
            vector<V> tempVec = theVec;
            theVec.clear();
            for (size_t i = 0; i < numRows; ++i)
            {
                if (!toDelete[i])
                {
                    theVec.push_back(tempVec[i]);
                }
            }
        }
    }


    void sortKeys()
    {
        numRows = mainArray.size();
        // Below: do not sort if there is nothing to do
        //    Also, if numRows = 0, numRows-1 is numeric_limits<size_t>::max().
        if (numRows <= 1)
        {
            return;
        }
        vector<size_t> tempPermutation(numRows);
        for (size_t i = 0; i < numRows; ++i)
        {
            tempPermutation[i] = i;
        }
        if (stable)
        {
            if (forward)
            {
                stable_sort(tempPermutation.begin(), tempPermutation.end(),
                            SortCompareForward(mainArray));
            }
            else
            {
                stable_sort(tempPermutation.begin(), tempPermutation.end(),
                            SortCompareReverse(mainArray));
            }
        }
        else
        {
            if (forward)
            {
                sort(tempPermutation.begin(), tempPermutation.end(),
                     SortCompareForward(mainArray));
            }
            else
            {
                sort(tempPermutation.begin(), tempPermutation.end(),
                     SortCompareReverse(mainArray));
            }
        }

        permutationVector.resize(numRows);
        for (size_t i = 0; i < numRows; ++i)
        {
            permutationVector[tempPermutation[i]] = i;
        }
        size_t keepFirstDelete = keepFirst ? 1 : 0;
        if (makeDistinct)
        {
            toDelete.resize(numRows, false);
            for (size_t i = 0; i < numRows - 1; ++i)
            {
                if (mainArray[tempPermutation[i]] ==
                        mainArray[tempPermutation[i+1]])
                {
                    toDelete[tempPermutation[i+keepFirstDelete]] = true;
                }
            }
            // Sort toDelete
            bool temp;
            tempPermutation = permutationVector;
            for (size_t i = 0; i < numRows; ++i)
            {
                while (tempPermutation[i] != i)
                {
                    // Cannot use swap for bit arrays in g++
                    temp = toDelete[i];
                    toDelete[i] = toDelete[tempPermutation[i]];
                    toDelete[tempPermutation[i]] = temp;
                    swap(tempPermutation[i], tempPermutation[tempPermutation[i]]);
                }
            }
        }
        tempPermutation.clear();
        apply(mainArray);
    }


    void setForward(bool wantForward)
    {
        forward = wantForward;
    }

    void setDistinct(bool distinctify)
    {
        makeDistinct = distinctify;
    }

    void setKeepFirst(bool theKeepFirst)
    {
        keepFirst = theKeepFirst;
    }

    void setStable(bool wantStable)
    {
        stable = wantStable;
    }

    SyncSort(vector<T>& inArray,
             bool wantForward = true,
             bool distinctify = false,
             bool theKeepFirst = true,
             bool wantStable = true) : mainArray(inArray)
    {
        setForward(wantForward);
        setDistinct(distinctify);
        setKeepFirst(theKeepFirst);
        setStable(wantStable);
        sortKeys();
    }


    void clear()
    {
        permutationVector.clear();
        if (makeDistinct)
        {
            toDelete.clear();
        }
    }


private:
    bool forward, makeDistinct, keepFirst, stable;
    size_t numRows;
    vector<T>& mainArray;
    vector<size_t> permutationVector;
    vector<bool> toDelete;

    // Below is the comparator we use for the array.
    class SortCompareForward
    {
    public:
        SortCompareForward(vector<T>& inputArray) : theArray(inputArray) { }
        bool operator()(size_t i, size_t j)
        {
            return (theArray[i] < theArray[j]);
        }
    private:
        vector<T>& theArray;
    };
    class SortCompareReverse
    {
    public:
        SortCompareReverse(vector<T>& inputArray) : theArray(inputArray) { }
        bool operator()(size_t i, size_t j)
        {
            return (!(theArray[i] < theArray[j]));
        }
    private:
        vector<T>& theArray;
    };
};

#endif

