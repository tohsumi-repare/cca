#ifndef MOLBIOLIB_SPARSEVECTOR_H
#define MOLBIOLIB_SPARSEVECTOR_H 1

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




/** \file SparseVector.hpp
 * Contains the #SparseVectorTraits and #SparseVector classese.
 */

/** The SparseVectorTraits class is used to define specific traits of
 * data types.  The below is for all things numeric.  Otherwise,
 * specializations are needed.
 */
template<typename T>
class SparseVectorTraits
{
public:
    static T nullValue()
    {
        return static_cast<T>(0);
    }
};
/** The SparseVectorTraits class is used to define specific traits of
 * data types.  The below is for strings.
 */
template<>
class SparseVectorTraits<string>
{
public:
    static string nullValue()
    {
        return "";
    }
};

/** The sparse vector is used in place of std::vector for sparse storage.
 * The only operations it has are initializations:
 * \code
 * SparseVector<long> x;
 * SparseVector<double> y(10);
 * SparseVector<string> z(20, "tester");
 * \endcode
 * <code>resize(size_t_value)</code>, <code>size()</code>, and access:
 * \code
 * y[2] = 25.7;  // The brackets are used when this is a lvalue or
 *               // access is of non-null elements that are present.
 * ...
 * double y2 = y.at(3);  // This is preferred when y is purely an
 *                             // rvalue.  Otherwise, if brackets are used and
 *                             // the index is not there, the entry will be
 *                             // created.
 * \endcode
 * Note that <code>resize</code> will erase any elements that have an index
 * larger than or equal to the new size.  Finally, one may erase an element by
 * referring to the index:
 * \code
 * y.erase(2);
 * \endcode
 */
template<typename T>
class SparseVector
{
public:

    SparseVector<T>& operator=(const SparseVector<T>& rhs)
    {
        if (this == &rhs)
        {
            return *this;
        }
        theData.clear();
        theData = rhs.theData;
        currSize = rhs.currSize;
        nullValue = rhs.nullValue;
        return *this;
    }


    // For strictly lvalue operations.
    T at(const size_t i)
    {
#ifdef PROG_DEBUG
        if (i >= currSize)
        {
            // Out of bounds error.
            assert(true == false);
        }
#endif
        if (theData.find(i) != theData.end())
        {
            return theData[i];
        }
        else
        {
            return nullValue;
        }
    }


    T& operator[](const size_t i)
    {
#ifdef PROG_DEBUG
        if (i >= currSize)
        {
            // Out of bounds error.
            assert(true == false);
        }
#endif
        // Create a new entry.
        if (theData.find(i) == theData.end())
        {
            theData[i] = nullValue;
        }
        return theData[i];
    }

    void erase(size_t i)
    {
        theData.erase(i);
    }

    void resize(size_t newSize)
    {
        if (newSize < currSize)
        {
            set<size_t> toDelete;
            for (typename unordered_map<size_t, T>::iterator i = theData.begin(); i != theData.end(); ++i)
                if ((i->first) >= newSize)
                {
                    toDelete.insert(i->first);
                }
            for (set<size_t>::iterator i = toDelete.begin(); i != toDelete.end(); ++i)
            {
                theData.erase(*i);
            }
        }
        currSize = newSize;
    }

    size_t size()
    {
        return currSize;
    }

    void push_back(T& newData)
    {
        ++currSize;
        theData[currSize - 1] = newData;
    }

    void clear()
    {
        theData.clear();
        currSize = 0;
    }


    SparseVector()
    {
        nullValue = SparseVectorTraits<T>::nullValue();
        currSize = 0;
    }
    SparseVector(size_t newSize)
    {
        nullValue = SparseVectorTraits<T>::nullValue();
        currSize = newSize;
    }
    SparseVector(size_t newSize, T newValue)
    {
        nullValue = SparseVectorTraits<T>::nullValue();
#ifdef DEBUG
        if (newValue != nullValue)
        {
            cerr << "Warning from SparseVector!  Initialization of a sparse vector is occuring with a non-null value.  This is contradictory to the usage of a sparse vector.  Proceeding with the execution." << endl;
        }
#endif
        currSize = newSize;
        for (size_t i = 0; i < newSize; ++i)
        {
            theData[i] = newValue;
        }
    }


    // Below is for when we use an unordered map.
    void setMaxLoadFactor(float maxLoadFactor)
    {
        theData.max_load_factor(maxLoadFactor);
    }


    // nullValue is left public in case someone wants to override the
    // null value to some other value than the default.
    T nullValue;
private:
    // map<size_t, T> theData;
    unordered_map<size_t, T> theData;
    size_t currSize;
};


#endif

