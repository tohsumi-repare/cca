#ifndef MOLBIOLIB_SYSTEM_H
#define MOLBIOLIB_SYSTEM_H 1

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


/** \file System.hpp
 * Call the system command.
 * Usage:
 * \code
 * string theCommand = "fill in command here";
 * int status = system(theCommand);
 * \endcode
 * This is used so that one can pass a string without doing a .c_str() call
 * and if PROG_DEBUG is defined, then this also prints out the system command.
 */
int system(string theCommand)
{
#ifdef PROG_DEBUG
    cerr << "System command: " << theCommand << endl;
#endif
    return system(theCommand.c_str());
}


#ifndef MOLBIOLIB_SYSTEM_STDOUT_MAX_BUFFER_SIZE
#define MOLBIOLIB_SYSTEM_STDOUT_MAX_BUFFER_SIZE 16777215
#endif
/** \file System.hpp
 * Call the system command and capture stdout output.
 * Modified and corrected from 
 * http://www.jeremymorgan.com/tutorials/c-programming/how-to-capture-the-output-of-a-linux-command-in-c/
 * Usage:
 * \code
 * string theCommand = "fill in command here";
 * string output = systemStdout(theCommand);
 * \endcode
 * This is used so that one can pass a string without doing a .c_str() call
 * and if PROG_DEBUG is defined, then this also prints out the system command.
 */
string systemStdout(string theCommand)
{
#ifdef PROG_DEBUG
    cerr << "System command: " << theCommand << endl;
#endif
    string result = "";
    FILE* stream;
    const int max_buffer = MOLBIOLIB_SYSTEM_STDOUT_MAX_BUFFER_SIZE;
    char* buffer = new char[MOLBIOLIB_SYSTEM_STDOUT_MAX_BUFFER_SIZE];
    theCommand.append(" 2>&1");

    stream = popen(theCommand.c_str(), "r");
    if (stream)
    {
        while (!feof(stream))
        {
            if (fgets(buffer, max_buffer, stream) != NULL)
            {
                result.append(buffer);
            }
        }
        pclose(stream);
    }
    return result;
}


#endif




