#ifndef MOLBIOLIB_FILEMODIFICATIONTIME_H
#define MOLBIOLIB_FILEMODIFICATIONTIME_H 1

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
#ifdef DEBUG
#include "src/Functions/SystemUtilities/FileSize.hpp"
#endif


/** \file FileModificationTime.hpp
 * Return a string of the time a file was modified in time_t.  This is the
 * number of seconds since midnight, 1st January 1970 GMT.  Smaller numbers
 * mean an older file.
 * Usage:
 * \code
 * string filename = "test.txt";
 * cout << fileModificationTime(filename) << endl;
 * \endcode
 */
time_t fileModificationTime(string filename)
{
#ifdef DEBUG
    if (fileSize(filename) == static_cast<ifstream::pos_type>(-1))
    {
        cerr << "Error in fileModificationTime.  " << filename << " does not exist.  Exiting." << endl;
        assert(true == false);
    }
#endif
    struct stat statbuf;
    stat(filename.c_str(), &statbuf);
    return statbuf.st_mtime;
}




/** \file FileModificationTime.hpp
 * Return true if filename1 is older than filename2.
 * Usage:
 * \code
 * string filename2 = "test.txt";
 * string filename1 = "test.txt.index";
 * // is filename1 older than filename2?
 * // Here, the if the index is older than the original file, regenerate index.
 * if (isFileOlderThan(filename1, filename2)) {
 *    // Regenerate index.
 * }
 * \endcode
 */
bool isFileOlderThan(string filename1, string filename2)
{
    // Let zeroTime = midnight, 1st January 1970 GMT
    if (fileModificationTime(filename1) <
            fileModificationTime(filename2))
    {
        return true;    // filename1 was created closer to zeroTime than filename2
    }
    else
    {
        return false;
    }
}

#endif

