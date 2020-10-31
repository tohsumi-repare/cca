#ifndef MOLBIOLIB_CONVERTVECTORSET_H
#define MOLBIOLIB_CONVERTVECTORSET_H 1

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


/** \file ConvertVectorSet.hpp
 * Convert between <code>vector</code> and <code>set</code>.
 */



/** \fn template<typename T> void convertVectorToSet(vector<T>& input, set<T>& output)
 * Convert from <code>vector</code> to <code>set</code>.  For example:
 * \code
 * vector<string> inVec;
 * // ... Populate inVec ...
 * set<string> outSet;
 * // The <string> below is optional, since the compiler can figure it out.
 * convertVectorToSet<string>(inVec, outSet);
 * \endcode
 */
template<typename T>
void convertVectorToSet(vector<T>& input, set<T>& output)
{
    output.clear();
    size_t n = input.size();
    for (size_t i = 0; i < n; ++i)
    {
        output.insert(input[i]);
    }
}


/** \fn template<typename T> set<T> convertVectorToSet(vector<T>& input)
 * Convert from <code>vector</code> to <code>set</code>.  For example:
 * \code
 * vector<string> inVec;
 * // ... Populate inVec ...
 * set<string> outSet;
 * // The <string> below is optional, since the compiler can figure it out.
 * outSet = convertVectorToSet<string>(inVec);
 * \endcode
 */
template<typename T>
set<T> convertVectorToSet(vector<T>& input)
{
    set<T> output;
    convertVectorToSet<T>(input, output);
    return output;
}


/** \fn template<typename T> void convertSetToVector(set<T>& input, vector<T>& output)
 * Convert from <code>set</code> to <code>vector</code>.  For example:
 * \code
 * set<string> inSet;
 * // ... Populate inSet ...
 * vector<string> outVec;
 * // The <string> below is optional, since the compiler can figure it out.
 * convertSetToVector<string>(inSet, outVec);
 * \endcode
 */
template<typename T>
void convertSetToVector(set<T>& input, vector<T>& output)
{
    output.clear();
    for (typename set<T>::iterator i = input.begin(); i != input.end(); ++i)
    {
        output.push_back(*i);
    }
}


/** \fn template<typename T> vector<T> convertSetToVector(set<T>& input)
 * Convert from <code>set</code> to <code>vector</code>.  For example:
 * \code
 * set<string> inSet;
 * // ... Populate inSet ...
 * vector<string> outVec;
 * // The <string> below is optional, since the compiler can figure it out.
 * outVec = convertSetToVector<string>(inSet);
 * \endcode
 */
template<typename T>
vector<T> convertSetToVector(set<T>& input)
{
    vector<T> output;
    convertSetToVector<T>(input, output);
    return output;
}


#endif


