#ifndef MOLBIOLIB_PRIMITIVETYPES_H
#define MOLBIOLIB_PRIMITIVETYPES_H 1

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


/** \file PrimitiveTypes.hpp
 * Define primitive types including nullptr.
 * All standard C++ includes are located in this file
 * so pathnames can be configured easily.
 */

#include <string>
// Need IO precision, so that setting precision in PrimitiveTypes makes sense.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

// Some files need this, even if DEBUG or PROG_DEBUG is not defined.
#include <cassert>

// Below is to support exceptions
#include <stdexcept>

// The below is needed for run-time type information (RTTI)
// where there is a typeid.
#include <typeinfo>

#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <cstdlib>

// Need below to define size_t
#include <cstddef>
// Need below to define numeric_limits<TYPE_HERE>::min(),::max(), etc.
#include <limits>

#include <math.h>
#include <cmath>

#include <vector>
#include <list>
#include <map>
#include <unordered_map>
#include <stack>
#include <queue>
#include <deque>
#include <tuple>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <iterator>
// Below defines #pair
#include <utility>

// For toupper()
#include <cctype>

// For isdigit()
#include <ctype.h>

// Below is for file modification times
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <ctime>

#include <sys/types.h>
#include <regex.h>
// G++'s regex does not work.  :-(  So, we use the POSIX version.
// Boost's doesn't compile under g++ 4.3.
// #include <regex>
// Below is so NULL is defined.  Typically, NULL should be defined in
// #include <cstddef> but just in case, we defining this manually.  We need
// this since regexec is a C routine (not C++) so must pass a C NULL construct
// instead of nullptr.
#ifndef NULL
// No defining of NULL to be ((void *)0) since in C++ not needed.
#define NULL 0
#endif


// Need to do below namespace, since may also need tr1 for some compilers.
using namespace std;


#if 0
/** nullptr class
 * Below is needed until switching to g++ 4.6 (where it is defined) or clang++.
 * The below is based on [Meyers96] S. Meyers. More
 * Effective C++, 2nd edition (Addison-Wesley, 1996).
 * The code is from the C++0x document SC22/WG21/N2431 = J16/07-0301
 * Change 1 to 0 on compilers with nullptr defined in cstddef.
 *
 * This is a const object...
 */
const class
{
public:
    // Below doxygen cannot parse these functions correctly so the conditional
    // macros will tell doxygen to ignore these methods.
    /// @cond nullptrPublic
    // convertible to any type of null non-member pointer...
    template<class T> operator T*() const
    {
        return 0;
    }
    // or any type of null member pointer...
    template<class C, class T> operator T C::*() const
    {
        return 0;
    }
    /// @endcond
private:
    /// @cond nullptrPrivate
    // whose address can't be taken
    void operator&() const;
    /// @endcond
} nullptr = {};   // and whose name is nullptr
#endif



/** \def MOLBIOLIB_IO_PRECISION
 * Define the number of digits of precision in file output for doubles.
 */
#define MOLBIOLIB_IO_PRECISION 15

/** \def HAVE_WC_AWK
 * Define below if system has wc and awk.  Speeds up the sequential
 * indexing of the ReadOnlyStringFile (and everything associated with it).
 * This has been tested only on GNU's wc and awk.  If not using GNU's tools,
 * one may need to edit src/Objects/ReadOnlyStringfile.hpp.
 * Otherwise, comment out below.  For example, one will need to comment out
 * the below for Mac OS X, since its version of wc does not have the -L option.
 */
#define HAVE_WC_AWK 1

/** \def RM_COMMAND
 * Define below for forced removal of a file (not directory).  This may
 * include wildcards.  Below rm is used in *NIX while del is used in DOS.
 */
#define RM_COMMAND "rm -f "
// #define RM_COMMAND "del "


#endif

