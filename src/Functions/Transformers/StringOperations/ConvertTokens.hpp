#ifndef MOLBIOLIB_CONVERTTOKENS_H
#define MOLBIOLIB_CONVERTTOKENS_H 1

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
#ifdef DEBUG
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/Transformers/TupleOperations/TokenizeTuple.hpp"
#include "src/Functions/Transformers/TupleOperations/GetTupleValues.hpp"
#endif



/** \file ConvertString.hpp
 * Convert <code>vector\<string\></code> to and from variadic data types.
 */


/** \fn template<typename... V> void convertTokensToArgs(vector<string>& tokens, V&... args)
 * Converts tokens to args.
 */
template<typename... V>
void convertTokensToArgs(vector<string>& tokens, V&... args)
{
    tuple<V...> tempTuple;
    tuplizeTokens(tokens, tempTuple);
    getTupleValues(tempTuple, args...);
}


/** \fn template<typename... V> void convertArgsToTokens(vector<string>& tokens, V&... args)
 * Converts variadic arguments to string tokens.
 * Note that due to the variadic nature, the tokens must be passed as the first
 * argument to this function, even though normally, it would be passed as the
 * second because it is the output variable.
 */
template<typename... V>
void convertArgsToTokens(vector<string>& tokens, V&... args)
{
    tuple<V...> tempTuple = make_tuple(args...);
    tokenizeTuple(tempTuple, tokens);
}

#endif

