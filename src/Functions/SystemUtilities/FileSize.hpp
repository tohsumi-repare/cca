#ifndef MOLBIOLIB_FILESIZE_H
#define MOLBIOLIB_FILESIZE_H 1

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




/** \file FileSize.hpp
 * Compute a file's size quickly, returning -1 if the file does not exist.
 * \fn ifstream::pos_type fileSize(string filename)
 * Computes the filesize given a filename.  Can be used to preallocate a
 * buffer size for the ReaderWriter's.
 * \section Usage Usage:
 * \code
 * ifstream::pos_type long fSize = fileSize("myGenome.fasta");
 * \endcode
 * Note that one may compare the output against -1 to check for file existence.
 */
ifstream::pos_type fileSize(string filename)
{

    ifstream fin(filename.c_str());
    if (!fin)
    {
        fin.close();
        // return numeric_limits<ifstream::pos_type>::max();
        // cannot be used since the above returns 0, which is not correct.
        return static_cast<ifstream::pos_type>(-1);
    }
    fin.clear();
    fin.seekg(0, ios_base::beg);
    ifstream::pos_type begin_pos = fin.tellg();
    fin.clear();
    fin.seekg(0, ios_base::end);
    ifstream::pos_type end_pos = fin.tellg();
    fin.close();
    return end_pos - begin_pos;
}

#endif
