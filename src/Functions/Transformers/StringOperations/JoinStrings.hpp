#ifndef MOLBIOLIB_JOINSTRINGS_H
#define MOLBIOLIB_JOINSTRINGS_H 1

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


// #Table (instead of #vector) was purposefully not used.  Most of the time,
// we are just splitting apart a line, so a #vector is the most natural.


/** \file JoinStrings.hpp
 * Join a vector of strings to form one string.
 * \fn void joinStrings(vector<string>& tokens, string delimiter, string& line)
 * This routine takes in a vector of strings and a character (which is
 * typically something like a
 * tab (i.e. &quot;\\t&quot;) or space &quot; &quot; and outputs
 * a strings which is the concatenation of strings in the vector separated by
 * the delimiter.
 * \section Usage Usage:
 * \code
 * vector<string> myTokens;
 * string outputLine = "";
 * // Fill myTokens.
 * joinStrings(myTokens, " | ", outputLine);
 * \endcode
 *
 */
void joinStrings(vector<string>& tokens, string delimiter, string& line)
{
    line = "";

    size_t tokensSize = tokens.size();
    if (tokensSize == 0)
    {
        return;
    }
    if (tokensSize > 1)
    {
        for (size_t j = 0; j < tokensSize - 1; ++j)
        {
            line += tokens[j] + delimiter;
        }
    }
    line += tokens[tokensSize - 1];
}



/** \file JoinStrings.hpp
 * Join a vector of strings to form one string.
 * \fn string joinStrings(vector<string>& tokens, string delimiter = "\t")
 * This routine takes in a vector of strings and a character (which is
 * by default a tab) and outputs
 * a strings which is the concatenation of strings in the vector separated by
 * the delimiter.
 * \section Usage Usage:
 * \code
 * vector<string> myTokens;
 * // Fill myTokens.
 * string outputLine = joinStrings(myTokens, " | ");
 * \endcode
 * where we separated by | instead of by tabs in the example.
 *
 */
string joinStrings(vector<string>& tokens, string delimiter = "\t")
{
    string line = "";
    joinStrings(tokens, delimiter, line);
    return line;
}


#endif

