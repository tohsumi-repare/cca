#ifndef MOLBIOLIB_FEATURESTABLEHUMANPROTEINSREADER_H
#define MOLBIOLIB_FEATURESTABLEHUMANPROTEINSREADER_H 1

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




/** \file FeaturesTableHumanProteinsReader.hpp
 * Read #Features from a UCSC blastHg18KG.txt file.
 * \fn void readFeaturesTableHumanProteins(SequenceFeatures& theFeatures, string proteinFile = "blastHg18KG.txt", bool sorted = false, SequenceFeatures::element_type startShift = 0, SequenceFeatures::element_type endShift = 0, size_t resLength = 0)
 * \code
 * SequenceFeatures newFeatures;
 * readFeaturesTableHumanProteins(newFeatures,
 *                               "blastHg18KG.txt", false, 0, 0, 0);
 * \endcode
 * The last seven parameters are optional.  By default, the protein file is
 * "blastHg18KG.txt".  Usually, the file is not sorted, so the default is
 * correct.  However, if set to <code>true</code>, then #Features can search for genes hit
 * at a particular #Location much faster.  The final 0 is the reserve length of the
 * string used to read in one line of rmsk.txt.  0 means not to reserve any
 * space and let the string class handle it.  Perhaps 500 is not a bad value
 * here.  The 0 and 0 are the start and stop shifts in coordinates, as in
 * #readFeaturesTable.
 *
 * Note that the query is the human genome.
 *
 */
void readFeaturesTableHumanProteins(SequenceFeatures& theFeatures,
                                    string proteinFile = "blastHg18KG.txt",
                                    bool sorted = false,
                                    SequenceFeatures::element_type startShift = 0,
                                    SequenceFeatures::element_type endShift = 0,
                                    size_t resLength = 0)
{


    theFeatures.setSorted(sorted);

    string line;
    if (resLength > 0)
    {
        line.reserve(resLength);
    }
    vector<string> tokens;

    ReadOnlyStringFile proteinTableFile(proteinFile, true, true, false, "readFeaturesTableHumanProtein.index", resLength);
    size_t currCount = 0;
    for (size_t i = 0; !proteinTableFile.fail(); ++i)
    {
        line = proteinTableFile[i];
        splitString(line, "\t", tokens);

        SequenceLocation newLoc;
        SequenceLocation::element_type theOrigin = newLoc.origin();
        newLoc.contig() = tokens[14];
        bool doSecondStrand = false;
        string theStrands = tokens[9];
        if (theStrands.size() == 1)
        {
            doSecondStrand = true;
            newLoc.strand() = "+";
        }
        else
        {
            newLoc.strand() = theStrands[1];    // Second strand is target's.
        }
        SequenceLocation::position_type newI;
        SequenceFeatures::element_type startPos = convertFromString<SequenceFeatures::element_type>(tokens[16]) + theOrigin;
        SequenceFeatures::element_type endPos = convertFromString<SequenceFeatures::element_type>(tokens[17]) + theOrigin;
        SequenceFeatures::element_type shiftedStart, shiftedEnd;
        if (newLoc.contig() == "+")
        {
            shiftedStart = startPos + startShift;
            shiftedEnd = endPos + endShift;
        }
        else
        {
            shiftedStart = startPos + endShift;
            shiftedEnd = endPos + startShift;
        }
        newI.assign(shiftedStart, shiftedEnd);
        newLoc.position() = newI;
        StringQualifiers newQualifiers;
        newQualifiers.push_back("matches", tokens[1]);
        newQualifiers.push_back("misMatches", tokens[2]);
        newQualifiers.push_back("repMatches", tokens[3]);
        newQualifiers.push_back("ncount", tokens[4]);
        newQualifiers.push_back("qNumInsert", tokens[5]);
        newQualifiers.push_back("qBaseInsert", tokens[6]);
        newQualifiers.push_back("tNumInsert", tokens[7]);
        newQualifiers.push_back("tBaseInsert", tokens[8]);
        newQualifiers.push_back("qStrand", convertToString<char>(theStrands[0]));  // query's strand.
        newQualifiers.push_back("qName", tokens[10]);   // Name of the
        // human protein
        // (accession number)
        newQualifiers.push_back("qSize", tokens[11]);
        newQualifiers.push_back("qStart", tokens[12]);
        newQualifiers.push_back("qEnd", tokens[13]);
        newQualifiers.push_back("tSize", tokens[15]);
        newQualifiers.push_back("blockCount", tokens[18]);
        newQualifiers.push_back("blockSizes", tokens[19]);
        newQualifiers.push_back("qStarts", tokens[20]);
        newQualifiers.push_back("tStarts", tokens[21]);
        newQualifiers.push_back("Gene Section", "Whole");
        // string currentIndex = tokens[10];
        // if (doSecondStrand)
        //    currentIndex = currentIndex + "+";
        string currentIndex = convertToString<size_t>(currCount);
        theFeatures.addFeature(currentIndex, newLoc, newQualifiers);
        ++currCount;
        if (doSecondStrand)
        {
            newLoc.strand() = "-";
            shiftedStart = startPos + endShift;
            shiftedEnd = endPos + startShift;
            newI.assign(shiftedStart, shiftedEnd);
            newLoc.position() = newI;
            // string currentIndex2 = tokens[10];
            // currentIndex2 = currentIndex2 + "-";
            string currentIndex2 = convertToString<size_t>(currCount);
            theFeatures.addFeature(currentIndex2, newLoc, newQualifiers);
            ++currCount;
        }
    }

    proteinTableFile.close();
}

#endif

