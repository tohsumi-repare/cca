#ifndef MOLBIOLIB_CSTDIO_H
#define MOLBIOLIB_CSTDIO_H 1

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
#include "src/Functions/SystemUtilities/System.hpp"


/** \file Cstdio.hpp
 * Give ability to pass strings to file commands in stdio.h.
 *
 * Usage:
 * \code
 * string filename1 = "fill in here";
 * string filename2 = "fill in here";
 * int status = remove(filename1);
 * status = rename(filename2, filename1);
 * \endcode
 * One may instead of a string in any of the above, also pass a const char *
 * to remove or rename.  However, in such a case, remove will not push
 * wildcard removals to system and will fail.
 * This is used so that one can pass a string without doing a .c_str() call.
 * If there is any wildcards in remove, then the system call
 *
 */


int remove(string filename)
{
    if ((filename.find(' ') != string::npos) ||
            ((filename.find('[') != string::npos) &&
             (filename.find(']') != string::npos)) ||
            (filename.find('*') != string::npos) ||
            (filename.find('?') != string::npos))
    {
        return system(RM_COMMAND + filename);
    }
    else
    {
        return remove(filename.c_str());
    }
}



int rename(string oldFilename, const char * newFilename)
{
    return rename(oldFilename.c_str(), newFilename);
}



int rename(const char * oldFilename, string newFilename)
{
    return rename(oldFilename, newFilename.c_str());
}


int rename(string oldFilename, string newFilename)
{
    return rename(oldFilename.c_str(), newFilename.c_str());
}


#endif

