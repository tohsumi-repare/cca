#ifndef MOLBIOLIB_TOKENIZETUPLE_H
#define MOLBIOLIB_TOKENIZETUPLE_H 1

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




template<size_t N, typename R, typename C>
class TokenTupleConvertor
{
public:
    static void convert(vector<string>& instrings, R& theRow)
    {
        istringstream is(instrings[N]);
        is >> get<N>(theRow);
        typedef typename tuple_element<N-1,R>::type Q;
        TokenTupleConvertor<N-1,R,Q>::convert(instrings, theRow);
    }
};
template<size_t N, typename R>
class TokenTupleConvertor<N,R,string>
{
public:
    static void convert(vector<string>& instrings, R& theRow)
    {
        get<N>(theRow) = instrings[N];
        typedef typename tuple_element<N-1,R>::type Q;
        TokenTupleConvertor<N-1,R,Q>::convert(instrings, theRow);
    }
};



template<typename R, typename C>
class TokenTupleConvertor<0,R,C>
{
public:
    static void convert(vector<string>& instrings, R& theRow)
    {
        istringstream is(instrings[0]);
        is >> get<0>(theRow);
    }
};

template<typename R>
class TokenTupleConvertor<0,R,string>
{
public:
    static void convert(vector<string>& instrings, R& theRow)
    {
        get<0>(theRow) = instrings[0];
    }
};


/** \file TokenizeTuple.hpp
 * Tokens <--> tuples.
 * \fn template<typename... T> void tuplizeTokens(vector<string>& instrings, tuple<T...>& theTuple)
 * Convert tokens to tuples.
 * \section Usage Usage:
 * \code
 * string line;  vector<string> tokens;
 * tuple<...put in items here...> theTuple;
 * #include "src/Functions/Transformers/StringOperations/StringSplit.hpp"
 * getline(fin, line);  // Below is if the input is a TSV
 * stringSplit(line, "\t", tokens);
 * tuplizeTokens(tokens, theTuple);
 * \endcode
 * Conversions are all via string streams.
 */
template<typename... T>
void tuplizeTokens(vector<string>& instrings, tuple<T...>& theTuple)
{
    typedef typename tuple_element< tuple_size< tuple<T...> >::value - 1,
            tuple<T...> >::type Q;
    TokenTupleConvertor< tuple_size<tuple<T...> >::value - 1,
                         tuple<T...>, Q>::convert(instrings, theTuple);
}






template<size_t N, typename R, typename C>
class TupleTokenConvertor
{
public:
    static void convert(vector<string>& outstrings, R& theRow)
    {
        ostringstream os;
        os << get<N>(theRow);
        outstrings[N] = os.str();
        typedef typename tuple_element<N-1,R>::type Q;
        TupleTokenConvertor<N-1,R,Q>::convert(outstrings, theRow);
    }
};
template<size_t N, typename R>
class TupleTokenConvertor<N,R,string>
{
public:
    static void convert(vector<string>& outstrings, R& theRow)
    {
        outstrings[N] = get<N>(theRow);
        typedef typename tuple_element<N-1,R>::type Q;
        TupleTokenConvertor<N-1,R,Q>::convert(outstrings, theRow);
    }
};

template<typename R, typename C>
class TupleTokenConvertor<0,R,C>
{
public:
    static void convert(vector<string>& outstrings, R& theRow)
    {
        ostringstream os;
        os << get<0>(theRow);
        outstrings[0] = os.str();
    }
};

template<typename R>
class TupleTokenConvertor<0,R,string>
{
public:
    static void convert(vector<string>& outstrings, R& theRow)
    {
        outstrings[0] = get<0>(theRow);
    }
};


/** \file TokenizeTuple.hpp
 * Tokens <--> tuples.
 * \fn template<typename... T> void tokenizeTuple(tuple<T...>& theTuple, vector<string>& outstrings)
 * Convert tuples to tokens.
 * \section Usage Usage:
 * \code
 * vector<string> tokens;
 * tuple<...put in items here...> theTuple;
 * ...
 * tokenizeTuple(theTuple, tokens);
 * \endcode
 * Conversions are all via string streams.
 */
template<typename... T>
void tokenizeTuple(tuple<T...>& theTuple, vector<string>& outstrings)
{
    typedef typename tuple_element< tuple_size< tuple<T...> >::value - 1,
            tuple<T...> >::type Q;
    outstrings.clear();
    outstrings.resize(tuple_size< tuple<T...> >::value);
    TupleTokenConvertor< tuple_size<tuple<T...> >::value - 1,
                         tuple<T...>, Q>::convert(outstrings, theTuple);
}

#endif

