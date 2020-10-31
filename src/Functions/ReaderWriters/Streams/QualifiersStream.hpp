#ifndef MOLBIOLIB_QUALIFIERSSTREAM_H
#define MOLBIOLIB_QUALIFIERSSTREAM_H 1

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


#include "src/Objects/Qualifiers.hpp"



/** \file QualifiersStream.hpp
 * This adds output streams to #Qualifiers.
 * \fn template<typename K, typename F> ostream& operator<<(ostream& out, Qualifiers<K,F>& theQualifiers)
 * This adds output streams to #Qualifiers.
 * \section Usage Usage:
 * \code
 * Qualifiers<string, string> theQualifiers;
 * cout << theQualifiers << endl;
 * \endcode
 * \section ProgNotes Programming Notes:
 * It is similar version of what's in Boost's example.
 * Note that there is no analog of an input stream since the definition of
 * the stream "[key1, [v1,v1a,v2,v3]]" would be ambiguous if the type of the
 * values allow commas in them, <i>e.g.</i> strings.
 */
template<typename K, typename F>
ostream& operator<<(ostream& out, Qualifiers<K,F>& theQualifiers)
{

    Table<K> theKeys = theQualifiers.getKeys();
    size_t numKeys = theKeys.size();
    out << "[";
    for (size_t i = 0; i < numKeys; ++i)
    {
        out << "[" << theKeys[i] << ",[";
        Table<F>& theValues = theQualifiers[theKeys[i]];
        size_t numValues = theValues.size();
        for (size_t j = 0; j < numValues; ++j)
        {
            out << theValues[j];
            if (j != numValues - 1)
            {
                out << ",";
            }
        }
        out << "]]";
        if (i != numKeys - 1)
        {
            out << ",";
        }
    }
    out << "]";

    return out;
}

#endif

