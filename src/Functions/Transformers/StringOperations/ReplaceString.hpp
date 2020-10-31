#ifndef MOLBIOLIB_REPLACESTRING_H
#define MOLBIOLIB_REPLACESTRING_H

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



/** \file ReplaceString.hpp
 * Replace all instances of a string.
 * \fn string replaceString(string theString, string replaceSpace = " ", string newSpace = "_")
 * By default, replaces all spaces by underscore.
 */
string replaceString(string theString, string replaceSpace = " ", string newSpace = "_")
{
    string newString = theString;
    size_t spaceSize = replaceSpace.size();
    size_t currPos = newString.find(replaceSpace);
    while (currPos != string::npos)
    {
        newString.replace(currPos, spaceSize, newSpace);
        currPos = newString.find(replaceSpace, currPos + 1);
    }
    return newString;
}



/** \file ReplaceString.hpp
 * Replace all instances of a string.
 * \fn string removeCharFromString(string inputStr, char theChar)
 * Removes all instances of theChar from a string inputStr.
 * This is a specialized case which is faster than using
 * <code> string stringChar(theChar); replaceString(inputStr, stringChar, "");
 * </code>.
 */
string removeCharFromString(string inputStr, char theChar)
{
    size_t numChars = 0;
    if (inputStr.size() == 0)
    {
        return "";
    }
    size_t strLen = inputStr.size();
    for (size_t i = 0; i < strLen; ++i)
    {
        if (inputStr[i] == theChar)
        {
            ++numChars;
        }
    }
    if (numChars == 0)
    {
        return inputStr;
    }
    string outputStr;
    outputStr.resize(strLen - numChars);
    size_t outPtr = 0;
    for (size_t i = 0; i < strLen; ++i)
    {
        if (inputStr[i] != theChar)
        {
            outputStr[outPtr] = inputStr[i];
            ++outPtr;
        }
    }
    return outputStr;
}


#endif

