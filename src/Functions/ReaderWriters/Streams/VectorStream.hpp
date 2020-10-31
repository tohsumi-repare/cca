#ifndef MOLBIOLIB_VECTORSTREAM_H
#define MOLBIOLIB_VECTORSTREAM_H 1

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
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Objects/SparseVector.hpp"

/** \file VectorStream.hpp
 * Adds stream operation to vectors.
 * Programming notes:
 * - convertToString and convertFromString functionality is
 *   replicated here, as vectorStreamsConvertToString and
 *   vectorStreamsConvertFromString, since convertToString and
 *   convertFromString may be used with a vector and thus give rise to a
 *   circular inclusion.
 * - These routines cannot handle a vector\< vector \<T\> \>
 *   (or more nested versions) type data structures as it will incur a
 *   self-defining problem.
 *
 */



/** \fn template<typename T> string vectorStreamsConvertToString(T theData)
 * Converts data to a <code>string</code> or #PrimarySequence.
 * See ConvertString.hpp for original function.
 */
template<typename T>
string vectorStreamsConvertToString(T theData)
{
    stringstream inputStream;
#ifdef DEBUG
    bool badConvert = (inputStream << theData).fail();
    if (badConvert)
    {
        cerr << "Error in vectorStreamsConvertToString!  Could not convert theData to a string.  Exiting since this only happens in patheological cases." << endl;
        assert(true == false);
    }
#else
    inputStream << theData;
#endif
    return inputStream.str();
}
/** \fn template<> string vectorStreamsConvertToString<bool>(bool theData)
 * Boolean specialization of #vectorStreamsConvertToString.
 * See ConvertString.hpp for original function.
 */
template<> string vectorStreamsConvertToString<bool>(bool theData)
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
/** \fn template<> string vectorStreamsConvertToString<string>(string theData)
 * String specialization of #vectorStreamsConvertToString.
 * See ConvertString.hpp for original function.
 */
template<> string vectorStreamsConvertToString<string>(string theData)
{
    return theData;
}

/** \fn template<typename T> T vectorStreamsConvertFromString(string theData)
 * Converts data to a <code>string</code> or #PrimarySequence.
 * See ConvertString.hpp for original function.
 */
template<typename T>
T vectorStreamsConvertFromString(string theData)
{
    T newData;
    stringstream outputStream(theData);
#ifdef DEBUG
    if (trimSpacesString(theData) == "")
    {
        cerr << "Warning in vectorStreamsConvertFromString.  Empty string passed.  Continuing execution." << endl;
    }
    bool badConvert = (outputStream >> newData).fail();
    if (badConvert)
    {
#ifndef PROG_DEBUG
        cerr << "Warning in vectorStreamsConvertFromString.  Could not convert " << theData << " to type " << typeid(newData).name() << ".  Continuing execution." << endl;
#else
        cerr << "Error in vectorStreamsConvertFromString.  Could not convert " << theData << " to type " << typeid(newData).name() << ".  Exiting." << endl;
        assert(true == false);
#endif
    }
#else
    outputStream >> newData;
#endif
    return newData;
}
/* \fn template<> bool vectorStreamsConvertFromString<bool>(string theData)
 * Boolean specialization of #vectorStreamsConvertFromString.
 * See ConvertString.hpp for original function.
 */
template<> bool vectorStreamsConvertFromString<bool>(string theData)
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
    cerr << "Error in vectorStreamsConvertFromString<bool>.  Unrecognized string.  Exiting.\n";
    assert(true == false);
}
/* \fn template<> ifstream::pos_type vectorStreamsConvertFromString<ifstream::pos_type>(string theData)
 * File position specialization of #vectorStreamsConvertFromString.
 * See ConvertString.hpp for original function.
 */
template<> ifstream::pos_type vectorStreamsConvertFromString<ifstream::pos_type>(string theData)
{
    // Need to do the below in two steps.  The conversion from string to
    // ifstream::pos_type seems to be ambiguous at least in g++ 4.3.1 and 4.4.4.
    size_t tempSize_t = vectorStreamsConvertFromString<size_t>(theData);
    return static_cast<ifstream::pos_type>(tempSize_t);
}

/* \fn template<> string vectorStreamsConvertFromString<string>(string theData)
 * String specialization of #vectorStreamsConvertFromString.
 * See ConvertString.hpp for original function.
 */
template<> string vectorStreamsConvertFromString<string>(string theData)
{
    return theData;
}



/** \file VectorStream.hpp
 * Adds stream operation to vectors.
 * \fn template<typename T> istream& operator>>(istream& in, vector<T>& theVector)
 * Adds an input stream operation to the vector class.  Input must be of
 * the form (&quot;a&quot; &quot;b&quot; &quot;c&quot; ... &quot;z&quot;)
 * where a through z are vector entries separated
 * by a space.
 * \todo Cannot handle the case of \&quot; in the token string.  Fix someday.
 */
template<typename T>
istream& operator>>(istream& in, vector<T>& theVector)
{
    string line = "";
    char c;
    in >> c;   // First character must be non-whitespace, i.e. (
#ifdef DEBUG
    assert (c == '(');
#endif
    while (c != ')')
    {
        in.get(c);   // Now, we read in the spaces too.
        if (c != ')')
        {
            line.push_back(c);
        }
    }

    string temp;
    vector<string> tokens;
    splitString(line, "\"", tokens);

    size_t numTokens = tokens.size();
    for (size_t i = 0; i < numTokens; ++i)
    {
        // Because there are double-quotes around everything,
        // the even-numbered entries are just spaces [between the quotes].
        if (i % 2 == 1)
        {
            theVector.push_back(vectorStreamsConvertFromString<T>(tokens[i]));
        }
    }

    return in;
}



/** \file VectorStream.hpp
 * Adds stream operation to vectors.
 * \fn template<typename T> ostream& operator<<(ostream& out, vector<T>& theVector)
 * Adds an output stream operation to the vector class.
 * Output is of the form
 * (&quot;a&quot; &quot;b&quot; &quot;c&quot; ... &quot;z&quot;)
 * where a through z are vector entries, separated by a space.
 * \todo Cannot handle the case of \&quot; in the token string.  Fix someday.
 */
template<typename T>
ostream& operator<<(ostream& out, vector<T>& theVector)
{
    out << '(';
    size_t numTokens = theVector.size();
    if (numTokens > 0)
    {
        for (size_t i = 0; i < numTokens - 1; ++i)
        {
            out << "\"" << vectorStreamsConvertToString<T>(theVector[i]) << "\" ";
        }
        out << "\"" << vectorStreamsConvertToString<T>(theVector[numTokens - 1]) << "\"";
    }
    out << ')';
    return out;
}


/** \file VectorStream.hpp
 * Adds generalized stream operation to vectors.
 * \fn template<typename Q, typename T> ostream& insertionOperator(ostream& out, vector<T>& theVector)
 * Adds an output stream operation to the vector class.
 * Output is of the form
 * a[delim]b[delim]c[delim]...[delim]z
 * where a through z are vector entries, separated by whatever delimiter is
 * used by Q.  No endl is produced.  An example usage is:
 * \code
 * vector<int> theData;
 * // fill theData
 * Ofstream ofp(filename);
 * // ...
 * ofp << "the data is = " << insertionOperator<CSVType, int>(theData)
 *     << "." << endl;
 * \endcode
 */
template<typename Q, typename T>
ostream& insertionOperator(ostream& out, vector<T>& theVector)
{
    Q joiner;
    size_t numTokens = theVector.size();

    if (numTokens > 0)
    {
        vector<string> tokens;
        for (size_t i = 0; i < numTokens-1; ++i)
        {
            out << vectorStreamsConvertToString<T>(theVector[i]) + joiner.delimiter();
        }
        out << vectorStreamsConvertToString<T>(theVector[numTokens-1]);
    }
    return out;
}


/** \file VectorStream.hpp
 * Adds stream operation to #SparseVector.
 * \fn template<typename T> istream& operator>>(istream& in, SparseVector<T>& theVector)
 * Adds an input stream operation to the #SparseVector class.  Input must be of
 * the form (&quot;a&quot; &quot;b&quot; &quot;c&quot; ... &quot;z&quot;)
 * where a through z are vector entries separated
 * by a space.
 * \todo Cannot handle the case of \&quot; in the token string.  Fix someday.
 */
template<typename T>
istream& operator>>(istream& in, SparseVector<T>& theVector)
{
    string line = "";
    char c;
    in >> c;   // First character must be non-whitespace, i.e. (
#ifdef DEBUG
    assert (c == '(');
#endif
    while (c != ')')
    {
        in.get(c);   // Now, we read in the spaces too.
        if (c != ')')
        {
            line.push_back(c);
        }
    }

    string temp;
    vector<string> tokens;
    splitString(line, "\"", tokens);

    T nullValue = SparseVectorTraits<T>::nullValue();
    T tempValue = nullValue;

    size_t numTokens = tokens.size();
    size_t currIndex = 0;
    for (size_t i = 0; i < numTokens; ++i)
    {
        // Because there are double-quotes around everything,
        // the even-numbered entries are just spaces [between the quotes].
        if (i % 2 == 1)
        {
            tempValue = vectorStreamsConvertFromString<T>(tokens[i]);
            if (tempValue != nullValue)
            {
                if (currIndex >= theVector.size())
                {
                    theVector.resize(currIndex+1);
                }
                theVector[currIndex] = tempValue;
            }
            ++currIndex;
        }
    }

    return in;
}


/** \file VectorStream.hpp
 * Adds stream operation to #SparseVector.
 * \fn template<typename T> ostream& operator<<(ostream& out, SparseVector<T>& theVector)
 * Adds an output stream operation to the vector class.
 * Output is of the form
 * (&quot;a&quot; &quot;b&quot; &quot;c&quot; ... &quot;z&quot;)
 * where a through z are vector entries, separated by a space.
 * \todo Cannot handle the case of \&quot; in the token string.  Fix someday.
 */
template<typename T>
ostream& operator<<(ostream& out, SparseVector<T>& theVector)
{
    out << '(';
    size_t numTokens = theVector.size();
    if (numTokens > 0)
    {
        // Using vector at(size_t) method since access is strictly lvalued.
        for (size_t i = 0; i < numTokens - 1; ++i)
        {
            out << "\"" << vectorStreamsConvertToString<T>(theVector.at(i)) << "\" ";
        }
        out << "\"" << vectorStreamsConvertToString<T>(theVector.at(numTokens - 1)) << "\"";
    }
    out << ')';
    return out;
}


/** \file VectorStream.hpp
 * Adds generalized stream operation to sparse vectors.
 * \fn template<typename Q, typename T> ostream& insertionOperator(ostream& out, SparseVector<T>& theVector)
 * Adds an output stream operation to the SparseVector class.
 * Output is of the form
 * a[delim]b[delim]c[delim]...[delim]z
 * where a through z are vector entries, separated by whatever delimiter is
 * used by Q.  No endl is produced.  An example usage is:
 * \code
 * SparseVector<int> theData;
 * // fill theData
 * Ofstream ofp(filename);
 * // ...
 * ofp << "the data is = " << insertionOperator<TSVType, int>(theData)
 *     << "." << endl;
 * \endcode
 */
template<typename Q, typename T>
ostream& insertionOperator(ostream& out, SparseVector<T>& theVector)
{
    Q joiner;
    size_t numTokens = theVector.size();

    if (numTokens > 0)
    {
        vector<string> tokens;
        for (size_t i = 0; i < numTokens-1; ++i)
        {
            out << vectorStreamsConvertToString<T>(theVector.at(i)) + joiner.delimiter();
        }
        out << vectorStreamsConvertToString<T>(theVector.at(numTokens-1));
    }
    return out;
}


#endif

