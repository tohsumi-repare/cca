#ifndef MOLBIOLIB_SPLITSTRING_H
#define MOLBIOLIB_SPLITSTRING_H 1

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
#include "src/Functions/Transformers/StringOperations/TrimSpacesString.hpp"



// #Table (instead of #vector) was purposefully not used.  Most of the time,
// we are just splitting apart a line, so a #vector is the most natural.


/** \file SplitString.hpp
 * Split and tokenize string.
 * \fn void splitString(string line, string delimiter, vector<string>& tokens, bool trimWhiteSpaces = false)
 * This routine takes in a string and a character (which is
 * typically something like a
 * tab (i.e. &quot;\\t&quot;) or space &quot; &quot; and outputs
 * a vector of strings which are the substrings before, between, and after
 * the delimiter.
 * \section Usage Usage:
 * \code
 * string inputLine = "Hello  there, world!";
 * vector<string> myTokens;
 * splitString(inputLine, " ", myTokens, true);
 * \endcode
 * where myTokens will have elements &quot;Hello&quot;,
 * &quot;there,&quot; and &quot;world!&quot;.  The
 * trailing optional parameter determines whether leading spaces and tabs
 * are trimmed off each entry.  By default, they are not trimmed.
 *
 * \section splitStringProgNotes Programmer's Notes
 * There are many versions of this routine out there.  Some of them (including
 * the last version of this) are flawed in that if a delimiter occured at the
 * end, a final null string is not pushed on the end of tokens.  This one is
 * not the fastest or most memory efficient, but one may put in a print
 * statement below to print out where all the delimiters were found.
 */
void splitString(string line, string delimiter, vector<string>& tokens, bool trimWhiteSpaces = false)
{

    tokens.clear();
    size_t delimSize = delimiter.size();

    vector<size_t> thePos;
    size_t foundPos = line.find(delimiter);
    while (foundPos != string::npos)
    {
        thePos.push_back(foundPos);
        foundPos = line.find(delimiter, foundPos+delimSize);
    }
    // *** Can put print statement here to print
    //     all locations of delimiters found ***
    size_t numTokensFound = thePos.size();
    size_t lastPos = 0;
    string tempToken;
    for (size_t i = 0; i < numTokensFound; ++i)
    {
        tempToken = line.substr(lastPos, thePos[i] - lastPos);
        if (trimWhiteSpaces)
        {
            tokens.push_back(trimSpacesString(tempToken));
        }
        else
        {
            tokens.push_back(tempToken);
        }
        lastPos = thePos[i] + delimSize;
    }
    tempToken = line.substr(lastPos, line.size() - lastPos);
    if (trimWhiteSpaces)
    {
        tokens.push_back(trimSpacesString(tempToken));
    }
    else
    {
        tokens.push_back(tempToken);
    }

}



/** \file SplitString.hpp
 * Split and tokenize string.
 * \fn size_t minSplitStringValueDelimIndex(vector<size_t>& values, size_t& index)
 * Returns both the minimum value of values as well as its index.
 * For use in splitString using vector of delimiters.
 */
size_t minSplitStringValueDelimIndex(vector<size_t>& values, size_t& index)
{
    size_t numValues = values.size();
    size_t minFoundPos = numeric_limits<size_t>::max();
    index = numeric_limits<size_t>::max();
    for (size_t i = 0; i < numValues; ++i)
    {
        if (values[i] < minFoundPos)
        {
            minFoundPos = values[i];
            index = i;
        }
    }
    return minFoundPos;
}

/** \file SplitString.hpp
 * Split and tokenize string.
 * \fn void splitString(string line, vector<string>& delimiter, vector<string>& tokens, bool trimWhiteSpaces = false, bool groupConsecutiveDelimters = false)
 * This routine is like the above, except it accepts a vector of delimiters.
 * It finds the position of the first delimiter it finds and splits there
 * unless <code>groupConsecutiveDelimters</code> is <code>true</code> in
 * which case it keeps looking until a non-delimiter or the end of the
 * string is found.
 * \todo Implement this some day.
 */
void splitString(string line, vector<string>& delimiter, vector<string>& tokens, bool trimWhiteSpaces = false, bool groupConsecutiveDelimters = false)
{

    tokens.clear();
    size_t numDelim = delimiter.size();
    vector<size_t> delimSize;
    for (size_t i = 0; i < numDelim; ++i)
    {
        delimSize.push_back(delimiter[i].size());
    }

    vector<size_t> thePos;
    vector<size_t> theDelim;
    vector<size_t> foundPos;
    for (size_t i = 0; i < numDelim; ++i)
    {
        foundPos.push_back(line.find(delimiter[i]));
    }
    size_t whichDelim = 0;
    size_t minFoundPos = minSplitStringValueDelimIndex(foundPos, whichDelim);
    while (minFoundPos != string::npos)
    {
        thePos.push_back(minFoundPos);
        theDelim.push_back(whichDelim);
        foundPos.clear();
        for (size_t i = 0; i < numDelim; ++i)
        {
            foundPos.push_back(line.find(delimiter[i],
                                         minFoundPos+delimSize[whichDelim]));
        }
        minFoundPos = minSplitStringValueDelimIndex(foundPos, whichDelim);
    }
    // *** Can put print statement here to print
    //     all locations of delimiters found ***
    size_t numTokensFound = thePos.size();
    size_t lastPos = 0;
    string tempToken = "", tempToken2 = "";
    for (size_t i = 0; i < numTokensFound; ++i)
    {
        tempToken = line.substr(lastPos, thePos[i] - lastPos);
        if (trimWhiteSpaces)
        {
            tempToken2 = trimSpacesString(tempToken);
            if (!groupConsecutiveDelimters || (tempToken2 != ""))
            {
                tokens.push_back(tempToken2);
            }
        }
        else
        {
            if (!groupConsecutiveDelimters || (tempToken != ""))
            {
                tokens.push_back(tempToken);
            }
        }
        lastPos = thePos[i] + delimSize[theDelim[i]];
    }
    tempToken = line.substr(lastPos, line.size() - lastPos);
    if (trimWhiteSpaces)
    {
        tempToken2 = trimSpacesString(tempToken);
        if (!groupConsecutiveDelimters || (tempToken2 != ""))
        {
            tokens.push_back(tempToken2);
        }
    }
    else
    {
        if (!groupConsecutiveDelimters || (tempToken != ""))
        {
            tokens.push_back(tempToken);
        }
    }
}


/** \file SplitString.hpp
 * Split and tokenize string.
 * \fn string beforeString(string org, string delimiter, bool findFirst = true, bool returnWhole = false)
 * This will return the part of the original string before the
 * delimiter.  An example usage is
 * \code
 * string temp = "gi|542225";
 * cout << beforeString(temp, "|", true, false);
 * \endcode
 * Note that the delimiter is <b>not</b> part of the output.  The
 * <code>true</code> (default) means to find the first instance of the
 * delimiter.  Otherwise, it finds the last instance.  The
 * <code>false</code> (default) means if the delimiter is not found,
 * an empty string is returned.  If <code>true</code>, then the original
 * word is returned.
 */
string beforeString(string old, string delimiter, bool findFirst = true, bool returnWhole = false)
{
    size_t pos;
    if (findFirst)
    {
        pos = old.find(delimiter);
    }
    else
    {
        pos = old.rfind(delimiter);
    }
    if (pos == string::npos)
    {
        if (returnWhole)
        {
            return old;
        }
        else
        {
            return "";
        }
    }
    return old.substr(0, pos);
}



/** \file SplitString.hpp
 * Split and tokenize string.
 * \fn string afterString(string org, string delimiter, bool findLast = true, bool returnWhole = false)
 * This will return the part of the original string after the
 * delimiter.  An example usage is
 * \code
 * string temp = "gi|542225";
 * cout << afterString(temp, "|", true, false);
 * \endcode
 * Note that the delimiter is <b>not</b> part of the output. The
 * <code>true</code> (default) means to find the last instance of the
 * delimiter.  Otherwise, it finds the first instance.  The
 * <code>false</code> (default) means if the delimiter is not found,
 * an empty string is returned.  If <code>true</code>, then the original
 * word is returned.
 */
string afterString(string old, string delimiter, bool findLast = true, bool returnWhole = false)
{
    size_t delimLen = delimiter.size();
    size_t pos;
    if (findLast)
    {
        pos = old.rfind(delimiter);
    }
    else
    {
        pos = old.find(delimiter);
    }
    if (pos == string::npos)
    {
        if (returnWhole)
        {
            return old;
        }
        else
        {
            return "";
        }
    }
    return old.substr(pos+delimLen);
}



/** \file SplitString.hpp
 * Split and tokenize string.
 * \fn void splitStringOnSpaces(string line, vector<string>& tokens)
 * Splits string on spaces using stringstream.
 * This is different from ssvString, since here any number of spaces may be
 * between entries.
 */
void splitStringOnSpaces(string line, vector<string>& tokens)
{
    tokens.clear();
    stringstream theStringStream(line);
    copy(istream_iterator<string>(theStringStream),
         istream_iterator<string>(),
         back_inserter< vector<string> >(tokens));
}


/** \file SplitString.hpp
 * Split and tokenize string.
 * \fn vector<string> splitStringOnSpaces(string line)
 * Splits string on spaces using stringstream.
 */
vector<string> splitStringOnSpaces(string line)
{
    vector<string> tokens;
    splitStringOnSpaces(line, tokens);
    return tokens;
}

#endif

