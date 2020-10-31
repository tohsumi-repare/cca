#ifndef MOLBIOLIB_FEATURESTABLEREFGENEREADER_H
#define MOLBIOLIB_FEATURESTABLEREFGENEREADER_H 1

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




/** \file FeaturesTableRefGeneReader.hpp
 * Read #Features from a UCSC refGene.txt file and an optional refLink.txt file.
 * \fn void readFeaturesTableRefGene(SequenceFeatures& theFeatures, string featureFile = "refGene.txt", string linkFile = "", bool sorted = false, bool collapseToStart = false, bool collapseToEnd = false, SequenceFeatures::element_type startShift = 0, SequenceFeatures::element_type endShift = 0, bool suppressExons = false, bool suppressIntrons = false, size_t resLength = 0)
 * \code
 * SequenceFeatures newFeatures;
 * readFeaturesTableRefGene(newFeatures,
 *                          "refGene.hg18.txt", "refLink.hg18.txt", false,
 *                          false, false, 0, 0, false, false, 0);
 * \endcode
 * The last seven parameters are optional.  By default, the feature file is
 * "refGene.txt".  The ref link file is missing, so is not read, but if it
 * is read, then additional qualifiers are available.  Usually, refGene.txt
 * is not sorted, so the default is correct.  However, if set
 * to <code>true</code>, then #Features can search for features hit at a particular
 * #Location much faster.  The last 0 is the reserve length of the string used to
 * read in one line of refGene.txt.  0 means not to reserve any space and let
 * the string class handle it.  Perhaps 1000 is not a bad value here.  The
 * 0 and 0 are the start and stop shifts as described in #readFeaturesTable.
 * These shifts are for the global feature positions only and not for the exon
 * coordinates.  The second and third <code>false</code> specifies
 * the global feature coordinates are not collapsed to the start and end
 * coordinates respectively.  Like #readFeaturesTableEnsemblStructures, this will
 * create one whole feature row as well as exon rows that can be differentiated by
 * the "Feature Section" qualifier ("Whole" or "Exon").  The third
 * <code>false</code> specifies that the featureration of exon lines will not be
 * suppressed.  The final <code>false</code> specifies that the featureration of
 * intron lines will not be suppressed.
 */
void readFeaturesTableRefGene(SequenceFeatures& theFeatures,
                              string featureFile = "refGene.txt",
                              string linkFile = "",
                              bool sorted = false,
                              bool collapseToStart = false,
                              bool collapseToEnd = false,
                              SequenceFeatures::element_type startShift = 0,
                              SequenceFeatures::element_type endShift = 0,
                              bool suppressExons = false,
                              bool suppressIntrons = false,
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

    unordered_map< string, tuple<string, string, string> > refLink;
    bool doRefLink = false;
    if (linkFile != "")
    {
        doRefLink = true;
        ReadOnlyStringFile refLinkTableFile(linkFile, true, true, false, "readFeaturesTableRefGene.index", resLength);
        for (size_t i = 0; !refLinkTableFile.fail(); ++i)
        {
            line = refLinkTableFile[i];
            splitString(line, "\t", tokens);
            tuple<string, string, string> newLink;
            get<0>(newLink) = tokens[0];  // refLink name
            get<1>(newLink) = tokens[1];  // product
            get<2>(newLink) = tokens[3];  // protAcc
            refLink[tokens[2]] = newLink;  // mrnaAcc
        }
        refLinkTableFile.close();
    }


    ReadOnlyStringFile refGeneTableFile(featureFile, true, true, false, "readFeaturesTableRefGene.index", resLength);
    size_t currentIndexI = 0;
    for (size_t i = 0; !refGeneTableFile.fail(); ++i, ++currentIndexI)
    {
        line = refGeneTableFile[i];
        splitString(line, "\t", tokens);

        SequenceLocation newLoc;
        SequenceLocation::element_type theOrigin = newLoc.origin();
        newLoc.contig() = tokens[2];
        newLoc.strand() = tokens[3];
        SequenceLocation::position_type newI;
        // RefGene files coordinates start from 0.
        SequenceFeatures::element_type startPos = convertFromString<SequenceFeatures::element_type>(tokens[4]) + theOrigin;
        SequenceFeatures::element_type endPos = convertFromString<SequenceFeatures::element_type>(tokens[5]) + theOrigin;
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
            startPos -= endShift;
            endPos -= startShift;
        }
        newI.assign(startPos, endPos);
        newLoc.position() = newI;
        StringQualifiers newQualifiers;
        newQualifiers.push_back("cdsStart", tokens[6]);
        newQualifiers.push_back("cdsEnd", tokens[7]);
        newQualifiers.push_back("exonCount", tokens[8]);
        newQualifiers.push_back("exonStarts", tokens[9]);
        newQualifiers.push_back("exonEnds", tokens[10]);
        newQualifiers.push_back("name2", tokens[12]);
        newQualifiers.push_back("cdsStartStat", tokens[13]);
        newQualifiers.push_back("cdsEndStat", tokens[14]);
        newQualifiers.push_back("exonFrames", tokens[15]);
        newQualifiers.push_back("mrnaAcc", tokens[1]);
        if (doRefLink)
        {
            tuple<string, string, string>& currLink = refLink[tokens[1]];
            newQualifiers.push_back("name", get<0>(currLink));
            newQualifiers.push_back("product", get<1>(currLink));
            newQualifiers.push_back("protAcc", get<2>(currLink));
        }
        string currentIndex = convertToString<size_t>(currentIndexI);
        StringQualifiers entireQualifiers = newQualifiers;
        StringQualifiers intronQualifiers = newQualifiers;
        entireQualifiers.push_back("Gene Section", "Whole");
        theFeatures.addFeature(currentIndex, newLoc, entireQualifiers);
        newQualifiers.push_back("Gene Section", "Exon");
        intronQualifiers.push_back("Gene Section", "Feature-wide Intron");
        size_t numExons = convertFromString<size_t>(tokens[8]);
        vector<string> exonStartCoords, exonEndCoords;
        splitString(tokens[9], ",", exonStartCoords);
        splitString(tokens[10], ",", exonEndCoords);
        for (size_t i = 0; i < numExons; ++i)
        {
            SequenceFeatures::element_type exonStartCoord = convertFromString<SequenceFeatures::element_type>(exonStartCoords[i]) + theOrigin;
            SequenceFeatures::element_type exonEndCoord = convertFromString<SequenceFeatures::element_type>(exonEndCoords[i]) + theOrigin;
            newI.assign(exonStartCoord, exonEndCoord);
            newLoc.position() = newI;
            ++currentIndexI;
            currentIndex = convertToString<size_t>(currentIndexI);
            if (!suppressExons)
            {
                theFeatures.addFeature(currentIndex, newLoc, newQualifiers);
            }
            if (!suppressIntrons && (numExons > 1) && (i < (numExons-1)))
            {
                SequenceFeatures::element_type exonNextStartCoord = convertFromString<SequenceFeatures::element_type>(exonStartCoords[i+1]);
                if ((exonNextStartCoord - exonEndCoord) > 2)
                {
                    newI.assign(exonEndCoord+1, exonNextStartCoord-1);
                }
                newLoc.position() = newI;
                ++currentIndexI;
                currentIndex = convertToString<size_t>(currentIndexI);
                theFeatures.addFeature(currentIndex, newLoc, intronQualifiers);
            }
        }
    }
    refGeneTableFile.close();

}

#endif

