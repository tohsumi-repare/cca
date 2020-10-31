#ifndef MOLBIOLIB_FEATURESTABLERMSKREADER_H
#define MOLBIOLIB_FEATURESTABLERMSKREADER_H 1

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


#include "src/Objects/Features.hpp"
#include "src/Objects/ReadOnlyStringFile.hpp"
#include "src/Objects/Location.hpp"
#include "src/Objects/Interval.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/ReaderWriters/Features/FeaturesTableRefGeneReader.hpp"





/** \file FeaturesTableRmskReader.hpp
 * Read #Features from a UCSC rmsk.txt file.
 * \fn void readFeaturesTableRmsk(SequenceFeatures& theFeatures, string rmskFile = "rmsk.txt", bool sorted = false, bool collapseToStart = false, bool collapseToEnd = false, SequenceFeatures::element_type startShift = 0, SequenceFeatures::element_type endShift = 0, size_t resLength = 0)
 * \code
 * SequenceFeatures newFeatures;
 * readFeaturesTableRmsk(newFeatures,
 *                       "rmsk.hg18.txt", false, false, false, 0, 0, 0);
 * \endcode
 * The last seven parameters are optional.  By default, the rmsk file is
 * "rmsk.txt".  Usually, rmsk.txt is not sorted, so the default is correct.
 * However, if set to <code>true</code>, then #Features can search for
 * genes hit
 * at a particular #Location much faster.  The final 0 is the reserve length of the
 * string used to read in one line of rmsk.txt.  0 means not to reserve any
 * space and let the string class handle it.  Perhaps 500 is not a bad value
 * here.  The 0 and 0 are the start and stop shifts in coordinates as
 * described in #readFeaturesTable.  The second and third <code>false</code> means
 * not to collapse the coordinate to the start and end coordinates repsectively.
 */
void readFeaturesTableRmsk(SequenceFeatures& theFeatures,
                           string rmskFile = "rmsk.txt",
                           bool sorted = false,
                           bool collapseToStart = false,
                           bool collapseToEnd = false,
                           SequenceFeatures::element_type startShift = 0,
                           SequenceFeatures::element_type endShift = 0,
                           size_t resLength = 0)
{

#ifdef DEBUG
    assert((collapseToStart && collapseToEnd) == false);
#endif

    theFeatures.setSorted(sorted);

    string line;
    if (resLength > 0)
    {
        line.reserve(resLength);
    }
    vector<string> tokens;

    ReadOnlyStringFile rmskTableFile(rmskFile, true, true, false, "readFeaturesTableRmsk.index", resLength);
    for (size_t i = 0; !rmskTableFile.fail(); ++i)
    {
        line = rmskTableFile[i];
        splitString(line, "\t", tokens);

        SequenceLocation newLoc;
        SequenceLocation::element_type theOrigin = newLoc.origin();

        newLoc.contig() = tokens[5];
        newLoc.strand() = tokens[9];
        SequenceLocation::position_type newI;
        SequenceFeatures::element_type startPos = convertFromString<SequenceFeatures::element_type>(tokens[6]) + theOrigin;
        SequenceFeatures::element_type endPos = convertFromString<SequenceFeatures::element_type>(tokens[7]) + theOrigin;
        if (newLoc.contig() == "+")
        {
            if (collapseToStart)
            {
                endPos = startPos;
            }
            if (collapseToEnd)
            {
                startPos = endPos;
            }
            startPos += startShift;
            endPos += endShift;
        }
        else
        {
            if (collapseToStart)
            {
                startPos = endPos;
            }
            if (collapseToEnd)
            {
                endPos = startPos;
            }
            startPos += endShift;
            endPos += startShift;
        }
        newI.assign(startPos, endPos);
        newLoc.position() = newI;
        StringQualifiers newQualifiers;
        newQualifiers.push_back("genoLeft", tokens[8]);
        newQualifiers.push_back("repName", tokens[10]);
        newQualifiers.push_back("repClass", tokens[11]);
        newQualifiers.push_back("repFamily", tokens[12]);
        newQualifiers.push_back("repStart", tokens[13]);
        newQualifiers.push_back("repEnd", tokens[14]);
        newQualifiers.push_back("repLeft", tokens[15]);
        newQualifiers.push_back("Gene Section", "Whole");

        string currentIndex = convertToString<size_t>(i);
        theFeatures.addFeature(currentIndex, newLoc, newQualifiers);
    }

    rmskTableFile.close();
}

#endif

