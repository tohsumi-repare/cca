#ifndef MOLBIOLIB_FILESTREAM_H
#define MOLBIOLIB_FILESTREAM_H 1

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



/** \file FileStream.hpp
 * Contains the #Ifstream and #Ofstream classes.
 */

/** Ifstream class
 * Usage: same as ifstream class, except can pass string instead of char *.
 *        open will also take in a string.
 */
class Ifstream : public ifstream
{
public:
    Ifstream() : ifstream() { }
    Ifstream( string filename, ios_base::openmode mode = ios_base::in ) : ifstream(filename.c_str(), mode) { }
    void open ( string filename, ios_base::openmode mode = ios_base::in )
    {
        ifstream::open(filename.c_str(), mode);
    }
};

/** Ofstream class
 * Usage: same as ofstream class, except can pass string instead of char *.
 *        Will also set the output precision of doubles to
 *        MOLBIOLIB_IO_PRECISION, defined in PrimitiveTypes.hpp.
 *        open will also take in a string.
 */
class Ofstream : public ofstream
{
public:
    Ofstream() : ofstream () { }
    Ofstream( const char * filename, ios_base::openmode mode = ios_base::out ) : ofstream(filename, mode)
    {
        this->precision(MOLBIOLIB_IO_PRECISION);
    }
    Ofstream( string filename, ios_base::openmode mode = ios_base::out ) : ofstream(filename.c_str(), mode)
    {
        this->precision(MOLBIOLIB_IO_PRECISION);
    }
    void open ( string filename, ios_base::openmode mode = ios_base::out )
    {
        ofstream::open(filename.c_str(), mode);
    }
};

#endif

