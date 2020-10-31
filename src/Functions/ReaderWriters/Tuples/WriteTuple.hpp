#ifndef MOLBIOLIB_WRITETUPLE_H
#define MOLBIOLIB_WRITETUPLE_H 1

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
#include "src/Functions/Transformers/TupleOperations/TokenizeTuple.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"


/** \file WriteTuple.hpp
 * Output a tuple as a delimited row.
 * \fn template<typename... T> void writeTuple(tuple<T...>& theTuple, ostream& out, string delimiter = "\t")
 * By default this writes to the ostream a TSV of a tuple.  However, the
 * delimiter may be changed by passing a delimiter string after the tuple.
 * Example usage is
 * \code
 * tuple<string, double> myTuple;
 * Ofstream myOut;
 * // ...
 * writeTuple(myTuple, myOut, ",");  // if want to do a CSV out instead.
 * \endcode
 */
template<typename... T>
void writeTuple(tuple<T...>& theTuple, ostream& out,
                string delimiter = "\t")
{
    vector<string> tokens;
    tokenizeTuple(theTuple, tokens);
    size_t numTokens = tokens.size();
    if (numTokens > 0)
    {
        for (size_t i = 0; i < numTokens - 1; ++i)
        {
            out << tokens[i] << delimiter;
        }
        out << tokens[numTokens - 1];
    }
}

#endif

