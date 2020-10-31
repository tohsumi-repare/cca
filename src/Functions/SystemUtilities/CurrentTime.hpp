#ifndef MOLBIOLIB_CURRENTTIME_H
#define MOLBIOLIB_CURRENTTIME_H 1

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

/** \file CurrentTime.hpp
 * Return a string of the current time.
 * Usage:
 * \code
 * cout << currentTime() << endl;
 * \endcode
 */
string currentTime()
{
    time_t calTime;
    time(&calTime);  // Gets the calendar time.
    tm *locTime = localtime(&calTime);  // Convert to local time.
    string theResult(asctime(locTime));
    if (theResult[theResult.size()-1] == '\n')
    {
        theResult.erase(theResult.size()-1);
    }
    return theResult;
}

#endif

