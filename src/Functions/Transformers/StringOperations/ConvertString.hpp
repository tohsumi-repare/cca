#ifndef MOLBIOLIB_CONVERTSTRING_H
#define MOLBIOLIB_CONVERTSTRING_H 1

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
#include "src/Functions/ReaderWriters/Streams/VectorStream.hpp"
#ifdef DEBUG
#include "src/Functions/Transformers/StringOperations/TrimSpacesString.hpp"
#endif



/** \file ConvertString.hpp
 * Convert <code>string</code>/#PrimarySequence to and from any data type.
 */


/** \fn template<typename T> string convertToString(T theData)
 * Converts data to a <code>string</code> or #PrimarySequence.  For example:
 * \code
 * double x = 3.14;
 * string pi;
 * pi = convertToString<double>(x);   // The <double> is not needed if it is
 *                                    // clear to the compiler the type.
 * \endcode
 */
template<typename T>
string convertToString(T theData)
{
    stringstream inputStream;
#ifdef DEBUG
    bool badConvert = (inputStream << theData).fail();
    if (badConvert)
    {
        cerr << "Error in convertToString!  Could not convert theData to a string.  Exiting since this only happens in patheological cases." << endl;
        assert(true == false);
    }
#else
    inputStream << theData;
#endif
    return inputStream.str();
}

/** \fn template<> string convertToString<long double>(long double theData)
 * Double specialization of #convertToString.  Needed to set precision.
 */
template<> string convertToString<long double>(long double theData)
{
    stringstream inputStream;
    inputStream.precision(MOLBIOLIB_IO_PRECISION);
#ifdef DEBUG
    bool badConvert = (inputStream << theData).fail();
    if (badConvert)
    {
        cerr << "Error in convertToString!  Could not convert theData to a string.  Exiting since this only happens in patheological cases." << endl;
        assert(true == false);
    }
#else
    inputStream << theData;
#endif
    return inputStream.str();
}

/** \fn template<> string convertToString<double>(double theData)
 * Double specialization of #convertToString.  Needed to set precision.
 */
template<> string convertToString<double>(double theData)
{
    stringstream inputStream;
    inputStream.precision(MOLBIOLIB_IO_PRECISION);
#ifdef DEBUG
    bool badConvert = (inputStream << theData).fail();
    if (badConvert)
    {
        cerr << "Error in convertToString!  Could not convert theData to a string.  Exiting since this only happens in patheological cases." << endl;
        assert(true == false);
    }
#else
    inputStream << theData;
#endif
    return inputStream.str();
}

/** \fn template<> string convertToString<float>(float theData)
 * Double specialization of #convertToString.  Needed to set precision.
 */
template<> string convertToString<float>(float theData)
{
    stringstream inputStream;
    inputStream.precision(MOLBIOLIB_IO_PRECISION);
#ifdef DEBUG
    bool badConvert = (inputStream << theData).fail();
    if (badConvert)
    {
        cerr << "Error in convertToString!  Could not convert theData to a string.  Exiting since this only happens in patheological cases." << endl;
        assert(true == false);
    }
#else
    inputStream << theData;
#endif
    return inputStream.str();
}

/** \fn template<> string convertToString<bool>(bool theData)
 * Boolean specialization of #convertToString.
 */
template<> string convertToString<bool>(bool theData)
{
    if (theData)
    {
        return "true";
    }
    else
    {
        return "false";
    }
}

/** \fn template<> string convertToString<string>(string theData)
 * String specialization of #convertToString.
 */
template<> string convertToString<string>(string theData)
{
    return theData;
}





/** \fn template<typename T> T convertFromString(string theData)
 * Converts data to a <code>string</code> or #PrimarySequence.  For example:
 * \code
 * PrimarySequence pi("3.14");
 * double x;
 * x = convertFromString<double>(pi);   // The <double> is needed since the
 *                                      // is not ambiguous due to implicit
 *                                      // casting.
 * \endcode
 */
template<typename T>
T convertFromString(string theData)
{
    T newData;
    stringstream outputStream(theData);
#ifdef DEBUG
    if (trimSpacesString(theData) == "")
    {
        cerr << "Warning in convertFromString.  Empty string passed.  Continuing execution." << endl;
    }
    bool badConvert = (outputStream >> newData).fail();
    if (badConvert)
    {
#ifndef PROG_DEBUG
        cerr << "Warning in convertFromString.  Could not convert " << theData << " to type " << typeid(newData).name() << ".  Continuing execution." << endl;
#else
        cerr << "Error in convertFromString.  Could not convert " << theData << " to type " << typeid(newData).name() << ".  Exiting." << endl;
        assert(true == false);
#endif
    }
#else
    outputStream >> newData;
#endif
    return newData;
}

/* \fn template<> string convertFromString<bool>(bool theData)
 * Boolean specialization of #convertToString.
 */
template<> bool convertFromString<bool>(string theData)
{
    if (theData == "true" ||
            theData == "True" ||
            theData == "TRUE" ||
            theData == "1")
    {
        return true;
    }
    else if (theData == "false" ||
             theData == "False" ||
             theData == "FALSE" ||
             theData == "0")
    {
        return false;
    }
    cerr << "Error in convertFromString<bool>.  Unrecognized string.  Exiting.\n";
    assert(true == false);
}

/* \fn template<> ifstream::pos_type convertFromString<string>(string theData)
 * File position specialization of #convertToString.
 */
template<> ifstream::pos_type convertFromString<ifstream::pos_type>(string theData)
{
    // Need to do the below in two steps.  The conversion from string to
    // ifstream::pos_type seems to be ambiguous at least in g++ 4.3.1 and 4.4.4.
    size_t tempSize_t = convertFromString<size_t>(theData);
    return static_cast<ifstream::pos_type>(tempSize_t);
}

/* \fn template<> string convertFromString<string>(string theData)
 * String specialization of #convertToString.
 */
template<> string convertFromString<string>(string theData)
{
    return theData;
}


#endif

