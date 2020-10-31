#ifndef MOLBIOLIB_TRIMSPACESSTRING_H
#define MOLBIOLIB_TRIMSPACESSTRING_H

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


/** \file TrimSpacesString.hpp
 * Trim white spaces (spaces and tabs) off ends of a string.
 */


/** \fn string trimSpacesString(string theString)
 * Trim white spaces (spaces and tabs and newlines) off ends of a string.
 * \code
 * String inTest = "  Tester ", outTest;
 * outTest = trimSpacesString(inTest);
 * \endcode
 * *** Warning!  This routine uses long for trimming off the end, so
 * technically, this routine can handle strings of maximum length
 * numeric_limits<long>::max(). ***
 * long is used since we have a -- that force long to be < 0, where as using
 * size_t would produce an unwanted wraparound.
 */
string trimSpacesString(string theString)
{
    // The below is so that there is not a case of wrap-around with size_t
    // being an unsigned integer.
    if (theString == "")
    {
        return "";
    }
    size_t start = 0;
    // end has to be of type long, since if end = 0 and did --end, then have
    // a wraparound for size_t of unsigned type.
    long end = static_cast<long>(theString.size()) - 1;
    while ((static_cast<long>(start) <= end) && (theString[start] == ' ' || theString[start] == '\t'))
    {
        ++start;
    }
    while ((end >= static_cast<long>(start)) && (theString[static_cast<size_t>(end)] == ' ' || theString[static_cast<size_t>(end)] == '\t' || theString[static_cast<size_t>(end)] == '\r' || theString[static_cast<size_t>(end)] == '\n'))
    {
        --end;
    }
    if (static_cast<long>(start) > end)
    {
        return "";
    }
    else
    {
        return theString.substr(start, static_cast<size_t>(end) - start + 1);
    }
}


#endif

