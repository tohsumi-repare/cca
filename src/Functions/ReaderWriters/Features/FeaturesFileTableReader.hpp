#ifndef MOLBIOLIB_FEATURESFILETABLEREADER_H
#define MOLBIOLIB_FEATURESFILETABLEREADER_H 1

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
#include "src/Objects/ReadOnlyDelimitedFile.hpp"
#include "src/Objects/Qualifiers.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"



/** \file FeaturesFileTableReader.hpp
 * Read #Features from a file (typically a tab-separated values (TSV) file).
 * \fn template<typename LocType> void readFileFeaturesTable(Features<LocType>& theFeatures, string filename, bool& hasStrand, string featuresExt = "featureStructure", typename LocType::position_type::element_type startShift = 0, typename LocType::position_type::element_type endShift = 0, bool collapseToStart = false, bool collapseToEnd = false, string defaultStrand = "+")
 * \code
 * Features<MyLocType> newFeatures;
 * readFileFeaturesTable(newFeatures, "myFeatures.tsv", "featureStructure");
 * \endcode
 * The optional third parameter is the extension to the features file (default
 * shown) which is used to get the structure of the features file.
 *
 * The featureStructure file is a tab separated values file consisting of a
 * keyword in the first column, a value in the second column, and for
 * qualifiers, the column number in the third column.  The keywords and their
 * defaults (if not provided) are:
 *
 *   - idPos = -1, no identification column if < 0.
 *   - contigColumn = 0, column for contigs.
 *   - strandColumn = -1, do strand column if < 0.
 *   - startColumn = 1, column for starting position of the feature.
 *   - endColumn = 2, column for ending position of the feature.
 *   - sorted = false, specify if the features file is sorted by chromosome and position.  Use sortFeaturesFile.cpp to sort, if needed.  Speeds up the search of a features search.
 *   - delimiter = "\t", split columns of feature file by this.
 *   - hasHeader = true, default says feature file has header.
 *   - trimWhiteSpaces = false, if true delete whitespaces around features.
 *   - fileIndex = 1, subtract this much off the positions in the start and end columns in the feature file.  Typically, to make a 1-based index file to the internal 0-based index, must subtract -1.  This is of type long, so can be negative.
 *   - senseStrand = +, the symbol(s) used for the sense strand in the feature file.
 *   - resLength = 0, if not 0, preallocate this amount for the longest line length in the feature.  Typically, not modified.
 *   - qualifier, after which there must be a qualifier's label and column number (0-index based).
 *
 * The above may be in appear in any order.  However, it is recommended that
 * for readability, the qualifiers are placed last in the
 * featuresStructure file.  An example of a file is given below for the UCSC
 * refGene.txt file:
 * \verbatim
idPos   1
contigColumn    2
strandColumn    3
startColumn     4
endColumn       5
sorted  false
hasHeader	false
qualifier       cdsStart        6
qualifier       cdsEnd  7
qualifier       exonCount       8
qualifier       exonStarts      9
qualifier       exonEnds        10
qualifier       name2   12
qualifier       cdsStartStat    13
qualifier       cdsEndStat      14
qualifier       exonFrames      15
qualifier       mrnaAcc 1 \endverbatim
 *
 * Notice that not all keywords need be specified in a featuresStructure file.
 */
template<typename LocType>
void readFileFeaturesTable(Features<LocType>& theFeatures, string filename,
                           bool& hasStrand,
                           string featuresExt = "featureStructure",
                           typename LocType::position_type::element_type startShift = 0,
                           typename LocType::position_type::element_type endShift = 0,
                           bool collapseToStart = false,
                           bool collapseToEnd = false,
                           string defaultStrand = "+")
{

#ifdef DEBUG
    assert((collapseToStart && collapseToEnd) == false);
#endif

    StringQualifiers theQualifiers;
    long idPos = -1;
    size_t contigColumn = 0;
    long strandColumn = -1;
    size_t startColumn = 1;
    size_t endColumn = 2;
    bool sorted = false;
    string senseStrand = "+";
    string delimiter = "\t";
    bool hasHeader = true;
    bool trimWhiteSpaces = false;
    long fileIndex = 1;
    size_t resLength = 0;

    string featuresFilename = filename + "." + featuresExt;
    ReadOnlyTSVFile fIfp(featuresFilename, true, true);
    vector<string> tokens;
    for (size_t i = 0; !fIfp.fail(); ++i)
    {
        tokens = fIfp[i];
        if (tokens[0] == "idPos")
        {
            idPos = convertFromString<long>(tokens[1]);
        }
        else if (tokens[0] == "contigColumn")
        {
            contigColumn = convertFromString<size_t>(tokens[1]);
        }
        else if (tokens[0] == "strandColumn")
        {
            strandColumn = convertFromString<long>(tokens[1]);
        }
        else if (tokens[0] == "startColumn")
        {
            startColumn = convertFromString<size_t>(tokens[1]);
        }
        else if (tokens[0] == "endColumn")
        {
            endColumn = convertFromString<size_t>(tokens[1]);
        }
        else if (tokens[0] == "sorted")
        {
            sorted = convertFromString<bool>(tokens[1]);
        }
        else if (tokens[0] == "delimiter")
        {
            delimiter = tokens[1];
        }
        else if (tokens[0] == "hasHeader")
        {
            hasHeader = convertFromString<bool>(tokens[1]);
        }
        else if (tokens[0] == "trimWhiteSpaces")
        {
            trimWhiteSpaces = convertFromString<bool>(tokens[1]);
        }
        else if (tokens[0] == "fileIndex")
        {
            fileIndex = convertFromString<long>(tokens[1]);
        }
        else if (tokens[0] == "senseStrand")
        {
            senseStrand = tokens[1];
        }
        else if (tokens[0] == "resLength")
        {
            resLength = convertFromString<size_t>(tokens[1]);
        }
        else if (tokens[0] == "qualifier")
        {
            // A qualifier line must have the keyword, the qualifier,
            // and the column number.
#ifdef DEBUG
            assert(tokens.size() >= 3);
#endif
            theQualifiers.push_back(tokens[1], tokens[2]);
        }
        else
        {
            cerr << "Warning in FeaturesTableReader!  Unknown token encountered in " << featuresFilename << ".  Ignoring it and continuing." << endl;
        }
    }
    fIfp.close();

    if (strandColumn < 0)
    {
        hasStrand = false;
    }
    else
    {
        hasStrand = true;
    }


    theFeatures.setSorted(sorted);

    ReadOnlyStringFile featuresTableFile(filename, true, true, hasHeader, "readFeaturesTableTSV.index", resLength);
    string line;
    if (resLength > 0)
    {
        line.reserve(resLength);
    }
    for (size_t i = 0; !featuresTableFile.fail(); ++i)
    {
        line = featuresTableFile[i];
        splitString(line, delimiter, tokens, trimWhiteSpaces);
        LocType newLoc;
        newLoc.contig() = tokens[contigColumn];
        if (strandColumn != -1)
        {
            newLoc.strand() = tokens[static_cast<size_t>(strandColumn)];
        }
        else
        {
            newLoc.strand() = defaultStrand;
        }
        typename Features<LocType>::position_type newI;
        typedef typename Features<LocType>::position_type::element_type ElemType;
        ElemType startPos = convertFromString<ElemType>(tokens[startColumn]) - fileIndex;
        ElemType endPos = convertFromString<ElemType>(tokens[endColumn]) - fileIndex;
        // Below: it is assumed endPos > startPos
        // So, on the antisense the endPos is actually the beginning.
        if ((strandColumn == -1) ||
                (newLoc.strand() == senseStrand))
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
        if (startPos > endPos)
        {
#ifdef DEBUG
            cerr << "Warning from readFeaturesTable.  startPos (" << startPos << ") is bigger than endPos (" << endPos << ").  Swapping them.\n";
#endif
            swap(startPos, endPos);
        }
        newI.assign(startPos, endPos);
        newLoc.position() = newI;
        StringQualifiers newQualifiers;

        // For each qualifier, get the appropriate column and
        // add it to the feature's description.
        if (!theQualifiers.empty())
        {
            for (StringQualifiers::iterator qualifierIterator = theQualifiers.begin();
                    qualifierIterator != theQualifiers.end(); ++qualifierIterator)
            {
                pair<string, string> currentQualifier = *qualifierIterator;
                int theColumn = convertFromString<int>(currentQualifier.second);
#ifdef DEBUG
                assert(static_cast<size_t>(theColumn) < tokens.size());
#endif
                newQualifiers.push_back(currentQualifier.first, tokens[theColumn]);
            }
        }

        if (idPos == -1)
        {
            theFeatures.addFeature(convertToString<size_t>(i), newLoc, newQualifiers);
        }
        else
        {
            theFeatures.addFeature(tokens[static_cast<size_t>(idPos)], newLoc, newQualifiers);
        }
    }
    featuresTableFile.close();
}





/**
 * The below is if we already know if the features has an strand field.
 */
template<typename LocType>
void readFileFeaturesTable(Features<LocType>& theFeatures, string filename,
                           string featuresExt = "featureStructure",
                           typename LocType::position_type::element_type startShift = 0,
                           typename LocType::position_type::element_type endShift = 0,
                           bool collapseToStart = false,
                           bool collapseToEnd = false,
                           string defaultStrand = "+")
{
    bool hasStrand = true;
    readFileFeaturesTable(theFeatures, filename, hasStrand,
                          featuresExt,
                          startShift, endShift,
                          collapseToStart, collapseToEnd,
                          defaultStrand);

}



#endif

