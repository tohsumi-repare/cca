#ifndef MOLBIOLIB_FEATURESBEDREADER_H
#define MOLBIOLIB_FEATURESBEDREADER_H 1

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


#include "src/Objects/FileStream.hpp"
#include "src/Functions/ReaderWriters/Features/FeaturesTableReader.hpp"



/** \file FeaturesBEDReader.hpp
 * Read #Features from a BED file.
 * \fn template<typename LocType> void readFeaturesTableBED(Features<LocType>& theFeatures, string filename, bool& hasStrand, bool sorted = false, size_t resLength = 0)
 * \code
 * Features<MyLocType> newFeatures;
 * bool hasStrand = true;
 * // Below will set hasStrand depending on what's in the BED file.
 * readFeaturesTableBED(newFeatures, "myFeatures.tsv", hasStrand, false, 0);
 * \endcode
 * The first <code>false</code> is to indicate that the coordinates in the
 * file are not sorted, i.e. features sorted first by contig, then sorted by
 * starting coordinate.  If <code>true</code>, then #Features can use a faster
 * method of returning the features associated with a particular location.
 * The 0 specifies that the input string length should not be
 * preallocated.  The latter should be specified for files where the line is
 * known to be long.
 *
 * Note that the below assumes that the "names" column of the BED may not
 * necessarily be unique, so creates its own ID number.
 */
template<typename LocType>
void readFeaturesTableBED(Features<LocType>& theFeatures, string filename, bool& hasStrand, bool sorted = false, size_t resLength = 0)
{
    StringQualifiers theQualifiers;

    // We have to do the below to get rid of the track lines, which do nothing.
    string filename2 = filename + ".noTrackLines";
    string theComm = "awk '$1 != \"track\"{print}' " + filename + " >& " + filename2;
    system(theComm);

    Ifstream ifp(filename2);
    string line = "";
    if (resLength > 0)
    {
        line.reserve(resLength);
    }
    while (!ifp.fail() && (line == ""))
    {
        getline(ifp, line);
    }
#ifdef DEBUG
    assert(line != "");
#endif
    vector<string> tokens;
    splitString(line, "\t", tokens);

    // Has name field.
    if (tokens.size() >= 4)
    {
        theQualifiers.push_back("name", "3");
    }
    if (tokens.size() >= 6)
    {
        hasStrand = true;
    }
    else
    {
        hasStrand = false;
    }
    if (hasStrand)
        readFeaturesTable(theFeatures, filename2, theQualifiers,
                          -1, 0, 5, 1, 2, sorted, "\t", false);
    else
        readFeaturesTable(theFeatures, filename2, theQualifiers,
                          -1, 0, -1, 1, 2, sorted, "\t", false);

    theComm = filename2 + " " + filename2 + "*.index";
    remove(theComm);
}

#endif

