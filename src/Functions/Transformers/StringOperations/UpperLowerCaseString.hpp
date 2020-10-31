#ifndef MOLBIOLIB_UPPERLOWERCASESTRING_H
#define MOLBIOLIB_UPPERLOWERCASESTRING_H 1

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


/** \file UpperLowerCaseString.hpp
 * Convert <code>string</code>/#PrimarySequence to and from upper/lower case.
 */


/** \fn string convertStringToUpperCase(string inputString)
 * Converts a <code>string</code> or #PrimarySequence to all upper case.
 * \code
 * string temp1 = "...";
 * string temp2 = convertStringToUpperCase(temp1);
 * \endcode
 */
string convertStringToUpperCase(string inputString)
{

    size_t n = inputString.size();
    string result = "";
    result.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        result[i] = toupper(inputString[i]);
    }
    return result;
}


/** \fn string convertStringToLowerCase(string inputString)
 * Converts a <code>string</code> or #PrimarySequence to all lower case.
 * \code
 * string temp1 = "...";
 * string temp2 = convertStringToLowerCase(temp1);
 * \endcode
 */
string convertStringToLowerCase(string inputString)
{

    size_t n = inputString.size();
    string result = "";
    result.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        result[i] = tolower(inputString[i]);
    }
    return result;
}

#endif

