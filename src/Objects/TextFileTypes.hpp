#ifndef MOLBIOLIB_TEXTFILETYPES_H
#define MOLBIOLIB_TEXTFILETYPES_H 1

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
#include "src/Functions/Transformers/StringOperations/CSVString.hpp"
#include "src/Functions/Transformers/StringOperations/SSVString.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/JoinStrings.hpp"


using namespace std;

// All of the below types must be of the form
// class NAME {
// public:
//    static string delimiter() {
//       // Return elimeter used.
//    }
//    static void splitLine(string line, vector<string>& tokens) {
//       // Code here to split a line into tokens.  Probably call
//       // splitString, csvString, or some other equivalent function.
//    }
//    static void joinTokens(vector<string>& tokens, string& line) {
//       // Code here to join tokens into a single line.  Probably call
//       // csvEntry for CSVs, otherwise perhaps a simple concatenation of
//       // strings.
//    }
// };
//
// Some of the types cannot use <code>joinStrings</code> since they require
// special processing of each element.

/** \file TextFileTypes.hpp
 * Defines specific text file types, e.g. CSV and TSV.
 * Defines the parsing behavior for comma-separated-values (CSV) files.
 */
class CSVType
{
public:
    static string delimiter()
    {
        return ",";
    }
    static void splitLine(string line, vector<string>& tokens)
    {
        bool trimWhiteSpaces = true;
        csvString(line, tokens, trimWhiteSpaces);
    }
    static void joinTokens(vector<string>& tokens, string& line)
    {
        line = "";
        if (tokens.size() == 0)
        {
            return;
        }
        if (tokens.size() > 1)
        {
            for (size_t j = 0; j < tokens.size() - 1; ++j)
            {
                line += csvEntry(tokens[j]) + ",";
            }
        }
        line += csvEntry(tokens[tokens.size() - 1]);
    }
};


/** \file TextFileTypes.hpp
 * Defines specific text file types, e.g. CSV and TSV.
 * Defines the parsing behavior for tab-separated-values (TSV) files.
 */
class TSVType
{
public:
    static string delimiter()
    {
        return "\t";
    }
    static void splitLine(string line, vector<string>& tokens)
    {
        bool trimWhiteSpaces = false;
        splitString(line, "\t", tokens, trimWhiteSpaces);
    }
    static void joinTokens(vector<string>& tokens, string& line)
    {
        joinStrings(tokens, "\t", line);
    }
};



/** \file TextFileTypes.hpp
 * Defines specific text file types, e.g. CSV and TSV.
 * Defines the parsing behavior for space-separated-values (SSV) files.
 */
class SSVType
{
public:
    static string delimiter()
    {
        return " ";
    }
    static void splitLine(string line, vector<string>& tokens)
    {
        ssvString(line, tokens);
    }
    static void joinTokens(vector<string>& tokens, string& line)
    {
        line = "";
        if (tokens.size() == 0)
        {
            return;
        }
        if (tokens.size() > 1)
        {
            for (size_t j = 0; j < tokens.size() - 1; ++j)
            {
                line += ssvEntry(tokens[j]) + " ";
            }
        }
        line += ssvEntry(tokens[tokens.size() - 1]);
    }
};



/** \file TextFileTypes.hpp
 * Defines specific text file types, e.g. CSV and TSV.
 * Defines the parsing behavior for spaces-separated-values files without
 * any quotes in the entries.  This is different from SSVType in that entries
 * may be separated by any number of spaces.
 */
class SpacesType
{
public:
    static string delimiter()
    {
        return " ";
    }
    static void splitLine(string line, vector<string>& tokens)
    {
        splitStringOnSpaces(line, tokens);
    }
    static void joinTokens(vector<string>& tokens, string& line)
    {
        // Note that this is not an exact inverse of splitLine since
        // the number of spaces between each entry is lost.
        joinStrings(tokens, " ", line);
    }
};




/** \file TextFileTypes.hpp
 * Defines specific text file types, e.g. CSV and TSV.
 * Defines the parsing behavior for single line files (StringTable).
 */
class StringType
{
public:
    static string delimiter()
    {
        return "";
    }
    static void splitLine(string line, vector<string>& tokens)
    {
        tokens.resize(1);
        tokens[0] = line;
    }
    static void joinTokens(vector<string>& tokens, string& line)
    {
        line = tokens[0];
    }
};




#endif

