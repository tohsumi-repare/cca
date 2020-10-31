#ifndef MOLBIOLIB_JOINTUPLES_H
#define MOLBIOLIB_JOINTUPLES_H 1

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




template<size_t N, size_t M, typename T, typename U>
class TuplesJoin
{
public:
    static void join(T& tupleIn, U& tupleOut)
    {
        get<N+M>(tupleOut) = get<N>(tupleIn);
        TuplesJoin<N-1,M,T,U>::join(tupleIn, tupleOut);
    }
};


template<size_t M, typename T, typename U>
class TuplesJoin<0,M,T,U>
{
public:
    static void join(T& tupleIn, U& tupleOut)
    {
        get<M>(tupleOut) = get<0>(tupleIn);
    }
};




/** \file JoinTuples.hpp
 * Join contents of two tuples.
 * \fn template<size_t M, typename... T, typename... U> void joinTuples(tuple<T...>& tupleIn, tuple<U...>& tupleOut)
 * \section Usage Usage:
 * \code
 * tuple<int, double> theTupleA;
 * tuple<string, char> theTupleB;
 * tuple<int, double, string, char> tupleOut;
 * joinTuples<2>(theTupleB, tupleOut);  // Put contents of tupleB into out
 *                                      // starting in the third column.
 *                                      // Does not affect other columns.
 * \endcode
 */
template<size_t M, typename... T, typename... U>
void joinTuples(tuple<T...>& tupleIn, tuple<U...>& tupleOut)
{
    TuplesJoin< tuple_size< tuple<T...> >::value - 1, M,
                tuple<T...>, tuple<U...> >::join(tupleIn, tupleOut);
}





/** \file JoinTuples.hpp
 * Join contents of two tuples.
 * \fn template<typename... T, typename... U> void joinTuples(tuple<T...>& tupleIn, tuple<U...>& tupleOut)
 * \section Usage Usage:
 * \code
 * tuple<int, double> theTupleA;
 * tuple<string, char> theTupleB;
 * tuple<int, double, string, char> tupleOut;
 * joinTuples(theTupleA, tupleOut);  // Put contents of tupleA into out
 *                                   // starting in the first column.
 *                                   // Does not affect other columns.
 * \endcode
 */
template<typename... T, typename... U>
void joinTuples(tuple<T...>& tupleIn, tuple<U...>& tupleOut)
{
    TuplesJoin< tuple_size< tuple<T...> >::value - 1, 0,
                tuple<T...>, tuple<U...> >::join(tupleIn, tupleOut);
}



/** \file JoinTuples.hpp
 * Join contents of two tuples.
 * \fn template<typename... T, typename... U, typename... V> void joinTuples(tuple<T...>& tupleA, tuple<U...>& tupleB, tuple<V...>& tupleOut)
 * \section Usage Usage:
 * \code
 * tuple<int, double> theTupleA;
 * tuple<string, char> theTupleB;
 * tuple<int, double, string, char> tupleOut;
 * joinTuples(theTupleA, theTupleB, tupleOut);  // Combine the above, where
 *                                              // tupleB's contents starts in
 *                                              // the column after tupleA's
 *                                              // last column in out.
 * \endcode
 */
template<typename... T, typename... U, typename... V>
void joinTuples(tuple<T...>& tupleA, tuple<U...>& tupleB,
                tuple<V...>& tupleOut)
{
    TuplesJoin< tuple_size< tuple<T...> >::value - 1, 0,
                tuple<T...>, tuple<V...> >::join(tupleA, tupleOut);
    TuplesJoin< tuple_size< tuple<U...> >::value - 1,
                tuple_size< tuple<T...> >::value,
                tuple<U...>, tuple<V...> >::join(tupleB, tupleOut);
}



#endif

