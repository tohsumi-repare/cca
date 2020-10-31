#ifndef MOLBIOLIB_CONVERTCSVSET_H
#define MOLBIOLIB_CONVERTCSVSET_H 1

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
#include "src/Functions/Transformers/VectorOperations/ConvertVectorSet.hpp"
#include "src/Objects/TextFileTypes.hpp"


/** \file ConvertCSVSet.hpp
 * Convert a CSV <code>string</code> to and from a <code>set</code>.
 */


/** \fn string convertSetToCSV(set<string>& inputSet)
 * Converts a string <code>set</code> to a CSV string.
 * \code
 * set<string> theSet;
 * // ... fill theSet.
 * string temp = convertSetToCSV(theSet);
 * \endcode
 */
string convertSetToCSV(set<string>& inputSet)
{
    string result = "";
    vector<string> tokens = convertSetToVector<string>(inputSet);
    CSVType::joinTokens(tokens, result);
    return result;
}


/** \fn set<string> convertCSVToSet(string& inputString);
 * Converts a CSV string to a string <code>set</code>.
 * \code
 * string theString = "...";
 * set<string> theSet = convertSetToCSV(theString);
 * \endcode
 */
set<string> convertCSVToSet(string& inputString)
{
    set<string> result;
    vector<string> tokens;
    CSVType::splitLine(inputString, tokens);
    convertVectorToSet<string>(tokens, result);
    return result;
}


#endif

