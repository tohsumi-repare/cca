#ifndef MOLBIOLIB_TUPLESTREAM_H
#define MOLBIOLIB_TUPLESTREAM_H 1

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



/** \file TupleStream.hpp
 * Adds stream operation to tuples.
 * \fn template<typename... T> istream& operator>>(istream& in, tuple<T...>& theTuple)
 * Adds an input stream operation to the #tuple class.  Input must be of
 * the form (&quot;a&quot; &quot;b&quot; &quot;c&quot; ... &quot;z&quot;)
 * where a through z are tuple entries separated
 * by a space.
 * \todo Cannot handle the case of \&quot; in the token string.  Fix someday.
 */
template<typename... T>
istream& operator>>(istream& in, tuple<T...>& theTuple)
{
    string line = "";
    char c;
    in >> c;   // First character must be non-whitespace, i.e. (
#ifdef DEBUG
    assert (c == '(');
#endif
    while (c != ')')
    {
        in.get(c);   // Now, we read in the spaces too.
        if (c != ')')
        {
            line.push_back(c);
        }
    }

    string temp;
    vector<string> pretokens, tokens;
    splitString(line, "\"", pretokens);

    size_t numTokens = pretokens.size();
    for (size_t i = 0; i < numTokens; ++i)
    {
        // Because there are double-quotes around everything,
        // the even-numbered entries are just spaces [between the quotes].
        if (i % 2 == 1)
        {
            tokens.push_back(pretokens[i]);
        }
    }
    tuplizeTokens(tokens, theTuple);

    return in;
}


/** \file TupleStream.hpp
 * Adds stream operation to #tuple's.
 * \fn template<typename... T> ostream& operator<<(ostream& out, tuple<T...>& theTuple)
 * Adds an output stream operation to the #tuple class.
 * Output is of the form
 * (&quot;a&quot; &quot;b&quot; &quot;c&quot; ... &quot;z&quot;)
 * where a through z are tuple entries, separated by a space.
 * \todo Cannot handle the case of \&quot; in the token string.  Fix someday.
 */
template<typename... T>
ostream& operator<<(ostream& out, tuple<T...>& theTuple)
{
    vector<string> tokens;
    tokenizeTuple(theTuple, tokens);
    out << '(';
    size_t numTokens = tokens.size();
    if (numTokens > 0)
    {
        for (size_t i = 0; i < numTokens - 1; ++i)
        {
            out << "\"" << tokens[i] << "\" ";
        }
        out << "\"" << tokens[numTokens - 1] << "\"";
    }
    out << ')';
    return out;
}

#endif

