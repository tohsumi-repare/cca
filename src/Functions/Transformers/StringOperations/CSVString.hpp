#ifndef MOLBIOLIB_CSVSTRING_H
#define MOLBIOLIB_CSVSTRING_H 1

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
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/Transformers/StringOperations/TrimSpacesString.hpp"


// In the below, we purposefully do not use a #Table.  Most of the time,
// we just need to split apart a string their components.

/** \file CSVString.hpp
 * Tokenize from a comma-separated-value (CSV) string.
 * \fn template<typename T> string csvEntry(T& entry)
 * This routine takes any type of data that can be streamed into a string
 * via <code>convertToString</code> and turns it into something that can
 * go into a CSV entry.  So, if there are leading or trailing whitespaces or
 * if there is a double-quote or comma within, then the entry will have
 * double quotes around them.  Double quotes within such an entry will be
 * turned into two double quotes, as per the standard.
 */
template<typename T>
string csvEntry(T& entry)
{
    string entryString = convertToString<T>(entry);
    size_t entrySize = entryString.size(),
           quotePos = entryString.find("\"", 0);
    if ((entrySize > 0 &&
            (entryString[0] == ' ' ||
             entryString[0] == '\t' ||
             entryString[entrySize-1] == ' ' ||
             entryString[entrySize-1] == '\t')) ||
            entryString.find(",", 0) != string::npos ||
            quotePos != string::npos)
    {
        if (quotePos != string::npos)
        {
            // Find all double quotes and replace them by two double quotes.
            string temp;
            size_t numQuotes = 0;
            for (size_t i = 0; i < entrySize; ++i)
                if (entryString[i] == '\"')
                {
                    ++numQuotes;
                }
            temp.resize(entrySize + numQuotes);
            for (size_t i = 0, tempPos = 0; i < entrySize; ++i, ++tempPos)
            {
                temp[tempPos] = entryString[i];
                if (entryString[i] == '\"')
                {
                    ++tempPos;
                    temp[tempPos] = '\"';
                }
            }
            entryString = temp;
        }
        return "\"" + entryString + "\"";
    }
    else
    {
        return entryString;
    }
}


/** \file CSVString.hpp
 * Tokenize from a comma-separated-value (CSV) string.
 * \fn string csvEntryToString(string entry, bool trimWhiteSpaces = false)
 * This takes a valid CSV entry which may have double quotes in them
 * and turns it into a "normal" string entry where the surrounding double
 * quotes are deleted and the two double quotes are replaced by a single quote.
 */
string csvEntryToString(string entry, bool trimWhiteSpaces = false)
{
    if (trimWhiteSpaces)
    {
        entry = trimSpacesString(entry);
    }
    if (entry.find('\"') != string::npos)
    {
        size_t firstPos = entry.find_first_of('\"'),
               lastPos = entry.find_last_of('\"'),
               entrySize = entry.size();
        // Preresult is the entry without the surrounding double-quotes.
        string preresult = entry.substr(0, firstPos);
        preresult += entry.substr(firstPos+1, lastPos - (firstPos+1));
        preresult += entry.substr(lastPos+1, entrySize - (lastPos+1));
        string result;
        result.reserve(preresult.size());
        result = "";
        // Turn all two double quotes into a single double quote.
        for (size_t i = 0; i < preresult.size(); ++i)
        {
            result += preresult[i];
            if (preresult[i] == '\"')
            {
                ++i;
            }
        }
        return result;
    }
    else
    {
        return entry;
    }
}



// Below simply determines if we are still in double quotes by counting the
// number of double quotes in tempEntry.  If the number of double quotes is
// odd, we are in the double quotes (return true).  This is on the assumption
// that the entry is well-formed.  For example, it is assumed all double-quotes
// in CSV entries between double quotes are of the form "" (two double quotes).
bool inCSVQuotes(string tempEntry)
{
    size_t entrySize = tempEntry.size(), count = 0;
    for (size_t i = 0; i < entrySize; ++i)
        if (tempEntry[i] == '\"')
        {
            ++count;
        }
    if (count % 2 == 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}


/** \file CSVString.hpp
 * Tokenize from a comma-separated-value (CSV) string.
 * \fn void csvString(string line, vector<string>& tokens, bool trimWhiteSpaces = false)
 * This routine takes in a string and outputs a  a vector of strings
 * that are the entries in the CSV string.  Note that a CSV has the
 * property that if any entry has a comma or double quote in it, the
 * entire entry is surrounded by double quotes.  Double quotes within
 * the entry are represented by "" (two double quotes).
 * \section Usage Usage:
 * \code
 * string inputLine = "Hello there, world!";
 * vector<string> myTokens;
 * csvString(inputLine, myTokens, false);
 * \endcode
 * where myTokens will have elements "Hello there" and " world!".  The optional
 * parameter, <code>trimWhiteSpaces</code>, allows one to specify if leading
 * and trailing spaces and tabs (before and after surrounding double quotes)
 * are to be trimmed off an entry.  The default is to trim.
 * \todo This is not quite the complete CSV parser.  A complete one would
 * also be able to handle multi-line input.  For example,
 * <pre>
 * Entry 1, "Entry 2 with a line
 * break in it"
 * </pre>
 */
void csvString(string line, vector<string>& tokens, bool trimWhiteSpaces = true)
{
    tokens.clear();

    size_t stringSize = line.size(), beginPos = 0, endPos;

    bool reached_end = false;

    while(beginPos < stringSize && !reached_end)
    {
        endPos = line.find(",", beginPos);
        if (endPos != string::npos)
        {
            if (endPos > beginPos)
            {
                // Below does not include the delimiter.
                while (endPos != string::npos &&
                        inCSVQuotes(line.substr(beginPos, endPos - beginPos)))
                {
                    endPos = line.find(",", endPos+1);
#ifdef DEBUG
                    if (endPos == string::npos &&
                            inCSVQuotes(line.substr(beginPos, stringSize - beginPos)))
                    {
                        cerr << "Error!  Invalid CSV input to csvString.  Exiting.\n";
                        assert(true == false);
                    }
#endif
                }
                if (endPos == string::npos)
                {
                    endPos = stringSize;
                }
                string tempS1 = line.substr(beginPos, endPos - beginPos);
                string tempS2 = csvEntryToString(tempS1, trimWhiteSpaces);
                tokens.push_back(tempS2);
                // tokens.push_back(csvEntryToString(line.substr(beginPos, endPos - beginPos), trimWhiteSpaces));
                beginPos = endPos + 1;
                if (endPos == stringSize - 1)
                {
                    reached_end = true;
                }
            }
            else
            {
                // Found comma at beginning.
                // I.e. a comma was followed by another.
                tokens.push_back("");
                beginPos++;
            }
        }
        else
        {
            // No delimiters found from beginPos to the end of the string.
            reached_end = true;
            tokens.push_back(csvEntryToString(line.substr(beginPos, stringSize - beginPos), trimWhiteSpaces));
        }
    }
}


/** \file CSVString.hpp
 * Tokenize from a comma-separated-value (CSV) string.
 * \fn void csvString(vector<string> line, vector<string>& tokens, bool trimWhiteSpaces = false)
 * This routine is like the csvString above except it takes a vector of strings.
 * This version can handle multi-line CSV entries.
 * \todo Implement this some day.
 */
void csvString(vector<string> lines, vector<string>& tokens, bool trimWhiteSpaces = true)
{
    assert(true == false);   // Not yet implemented.
}

#endif

