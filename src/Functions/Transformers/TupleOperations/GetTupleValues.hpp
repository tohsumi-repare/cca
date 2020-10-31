#ifndef MOLBIOLIB_GETTUPLEVALUES_H
#define MOLBIOLIB_GETTUPLEVALUES_H 1

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

// Below is an internal function, no to be used by any other
// function or object other than getTupleValues
template<size_t N, typename... V>
void getTupleValuesMolBioLibInternal(tuple<V...>& values)
{
    // Stop when there are no more arguments.
}

// Below is an internal function, no to be used by any other
// function or object other than getTupleValues
template<size_t N, typename... V, typename T, typename... Args>
void getTupleValuesMolBioLibInternal(tuple<V...>& values,
                                     T& currArg, Args&... args)
{
    currArg = get<tuple_size< tuple<V...> >::value - N>(values);
    getTupleValuesMolBioLibInternal<N-1, V...>(values, args...);
}


/** \file GetTupleValues.hpp
 * Get values of tuples in a variadic set of arguments
 * \fn template<typename... V> void getTupleValues(tuple<V...>& values, V&... args)
 * \section Usage Usage:
 * \code
 * tuple<string, string, double, string, int> testTuple;
 * // Fill testTuple...
 * string s1 = "", s2 = "", s3 = "";
 * double d1 = 0.0;
 * int i1 = 0;
 * // It is important that the order of variables below match
 * // the order of the tuple's declaraction's type .
 * getTupleValues(testTuple, s1, s2, d1, s3, i1);
 * \endcode
 *
 * Note that C++11 already has the inverse of this function which makes a
 * tuple from a set of arguments, called <code>make_tuple</code>:
 * \code
 * testTuple = make_tuple("xt1", "xt2", 3.14, "xt3", 6);
 * // or one can use variables...
 * testTuple = make_tuple(s1, s2, d1, s3, i1);
 * \endcode
 */

template<typename... V>
void getTupleValues(tuple<V...>& values, V&... args)
{
    getTupleValuesMolBioLibInternal<tuple_size< tuple<V...> >::value, V...>(values, args...);
}


#endif

