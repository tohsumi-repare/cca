#ifndef MOLBIOLIB_FEATURESTABLEENSEMBLQUALIFIERSREADER_H
#define MOLBIOLIB_FEATURESTABLEENSEMBLQUALIFIERSREADER_H 1

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



/** \file FeaturesTableEnsemblQualifiersReader.hpp
 * Read #Features from a mart_export.txt file from Ensembl's Biomart.
 * \fn void readFeaturesTableEnsemblQualifiers(SequenceFeatures& theFeatures, string featureFile = "mart_export.txt", bool sorted = false, bool collapseToStart = false, bool collapseToEnd = false, SequenceFeatures::element_type startShift = 0, SequenceFeatures::element_type endShift = 0, size_t resLength = 0)
 * You must
 * select <em>all</em> of the Feature attributes in order - first the left column,
 * then the right.  Example usage is:
 * \code
 * SequenceFeatures newFeatures;
 * readFeaturesTableEnsemblQualifiers(newFeatures,
 *                                    "mart_export.hg18.txt", false,
 *                                    false, false, 0, 0, 0);
 * \endcode
 * The last seven parameters are optional.  By default, the feature file is
 * "mart_export.txt".  mart_export.txt is not sorted, so the default is
 * correct.  However, if set to <code>true</code>, then #Features can search for
 * features hit at a particular #Location much faster.  The final 0 is the reserve
 * length of the string used to read in one line of mart_export.txt.  0 means
 * not to reserve any space and let the string class handle it.  Perhaps 2000 is
 * not a bad value here.  One may add a shift to the
 * global start and stop coordinates (0 and 0), like #readFeaturesTable.  The
 * last two <code>false</code> means not
 * to collapse the coordinates to the start and end coordinates respectively.
 */
void readFeaturesTableEnsemblQualifiers(SequenceFeatures& theFeatures,
                                        string featureFile = "mart_export.txt",
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

    // Ensembl has a header.
    ReadOnlyStringFile ensemblTableFile(featureFile, true, true, true, "readFeaturesTableEnsemblQualifiers.index", resLength);
    for (size_t i = 0; !ensemblTableFile.fail(); ++i)
    {
        line = ensemblTableFile[i];
        splitString(line, "\t", tokens);

        SequenceLocation newLoc;
        SequenceLocation::element_type theOrigin = newLoc.origin();
        newLoc.contig() = tokens[5];
        if (tokens[8] == "1")
        {
            newLoc.strand() = "+";
        }
        else if (tokens[8] == "-1")
        {
            newLoc.strand() = "-";
        }
        SequenceLocation::position_type newI;
        SequenceFeatures::element_type startPos = convertFromString<SequenceFeatures::element_type>(tokens[6]) - 1 + theOrigin;
        SequenceFeatures::element_type endPos = convertFromString<SequenceFeatures::element_type>(tokens[7]) - 1 + theOrigin;
        if (newLoc.strand() == "+")
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
            startPos -= endShift;
            endPos -= startShift;
        }
        newI.assign(startPos, endPos);
        newLoc.position() = newI;
        StringQualifiers newQualifiers;
        newQualifiers.push_back("Ensembl Feature ID", tokens[0]);
        newQualifiers.push_back("Ensembl Transcript ID", tokens[1]);
        newQualifiers.push_back("Ensembl Protein ID", tokens[2]);
        newQualifiers.push_back("Canonical transcript stable ID(s)", tokens[3]);
        newQualifiers.push_back("Description", tokens[4]);
        newQualifiers.push_back("Band", tokens[9]);
        newQualifiers.push_back("Transcript Start (bp)", tokens[10]);
        newQualifiers.push_back("Transcript End (bp)", tokens[11]);
        newQualifiers.push_back("Associated Feature Name", tokens[12]);
        newQualifiers.push_back("Associated Transcript Name", tokens[13]);
        newQualifiers.push_back("Associated Feature DB", tokens[14]);
        newQualifiers.push_back("Associated Transcript DB", tokens[15]);
        newQualifiers.push_back("Transcript count", tokens[16]);
        newQualifiers.push_back("% GC content", tokens[17]);
        newQualifiers.push_back("Feature Biotype", tokens[18]);
        newQualifiers.push_back("Transcript Biotype", tokens[19]);
        newQualifiers.push_back("Source", tokens[20]);
        newQualifiers.push_back("Status (gene)", tokens[21]);
        newQualifiers.push_back("Status (transcript)", tokens[22]);
        newQualifiers.push_back("Gene Section", "Whole");
        string currentIndex = convertToString<size_t>(i);
        theFeatures.addFeature(currentIndex, newLoc, newQualifiers);
    }

    ensemblTableFile.close();
}

#endif

