#ifndef MOLBIOLIB_FEATURESTABLEREADER_H
#define MOLBIOLIB_FEATURESTABLEREADER_H 1

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
#include "src/Objects/Qualifiers.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"



/** \file FeaturesTableReader.hpp
 * Read #Features from a file (typically a tab-separated values (TSV) file).
 * \fn template<typename LocType> void readFeaturesTable(Features<LocType>& theFeatures, string filename, StringQualifiers& theQualifiers, long idPos = -1, size_t contigColumn = 0, long strandColumn = -1, size_t startColumn = 1, size_t endColumn = 2, bool sorted = false, string delimiter = "\t", bool haveHeader = true, bool trimWhiteSpaces = false, typename LocType::position_type::element_type startShift = 0, typename LocType::position_type::element_type endShift = 0, bool collapseToStart = false, bool collapseToEnd = false, size_t resLength = 0, string defaultStrand = "+", long fileIndex = 0, string senseStrand = "+" )
 * \code
 * Features<MyLocType> newFeatures;
 * StringQualifiers theQualifiers;
 * // In the below, the column numbers (0-index based) of the qualifiers
 * // must be in double quotes.
 * theQualifiers.push_back("classification", "8");
 * theQualifiers.push_back("protein", "12");
 * readFeaturesTable(newFeatures, "myFeatures.tsv",
 *                   theQualifiers,
 *                   -1, 0, -1, 1, 2, false, "\t", true, false,
 *                   0, 0, 0, false, false, "+", 0, "+");
 * \endcode
 * All of the optional parameters are on the third and fourth lines and are
 * shown with their default values.  One may optionally not pass
 * theQualifiers, meaning there are no additional Qualifiers.
 * The 0 is the column (0-index) of the unique feature identifier where if -1
 * then use a simple row count as the feature identifier,
 * 0 is the contig column, the second -1 is the strand column,
 * 1 is the start position column, and 2 is the end position column.
 * <em>It is assumed the start position is <= the end position.</em>
 * If the strand column is specified to be -1, then the strand is not assumed
 * to be present in the file and is not read and all features have the strand
 * denoted by the default strand string (the first parameter, "+").  The
 * string used to denote the sense strand is give by the second "+" parameter.
 * The first <code>false</code> is to indicate that the coordinates in the
 * file are not sorted, i.e. features sorted first by contig, then sorted by
 * starting coordinate.  If <code>true</code>, then #Features can use a faster
 * method of returning the features associated with a particular location.
 * The tab is the delimiter between columns, the <code>true</code>
 * specifies the file has a header line, the second <code>false</code>
 * specifies that the white spaces around entries are not to be trimmed.
 * One may specify a global addition to the start and end boundary position of
 * the feature by specifying something other than 0 and 0 (last two 0's).  On the
 * antisense, it is a subtraction unless the strand is ignored.  Note
 * that start means the 5' end of the feature (so right side for features on the
 * antisense stand).  One can specify a negative (for subtracting) or positive
 * number.  The last two falses means do not collapse the coordiantes to the
 * start and stop coordinates respectively.  The point of using these is so
 * that one can use that in conjunction to the shifts to get the promoter or
 * UTR regions.  The second to last 0 specifies that the input string length should not be
 * preallocated.  The latter should be specified for files where the line is
 * known to be long.  The second to last parameter is the index base of the
 * file, with the default as 0.
 *
 * <em>Note that we do not specify the qualifier "Feature Section" in this reader.
 * However, in other readers an entire feature is denoted by "Whole", so perhaps
 * one may wish to include that in the qualifiers being passed.</em>
 */
template<typename LocType>
void readFeaturesTable(Features<LocType>& theFeatures, string filename,
                       StringQualifiers& theQualifiers,
                       long idPos = -1,
                       size_t contigColumn = 0,
                       long strandColumn = -1,
                       size_t startColumn = 1,
                       size_t endColumn = 2,
                       bool sorted = false,
                       string delimiter = "\t",
                       bool haveHeader = true,
                       bool trimWhiteSpaces = false,
                       typename LocType::position_type::element_type startShift = 0,
                       typename LocType::position_type::element_type endShift = 0,
                       bool collapseToStart = false,
                       bool collapseToEnd = false,
                       size_t resLength = 0,
                       string defaultStrand = "+",
                       long fileIndex = 0,
                       string senseStrand = "+" )
{

#ifdef DEBUG
    assert((collapseToStart && collapseToEnd) == false);
#endif

    theFeatures.setSorted(sorted);

    ReadOnlyStringFile featuresTableFile(filename, true, true, haveHeader, "readFeaturesTableTSV.index", resLength);
    string line;
    if (resLength > 0)
    {
        line.reserve(resLength);
    }
    vector<string> tokens;
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
                (newLoc.contig() == senseStrand))
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
                if (static_cast<size_t>(theColumn) >= tokens.size())
                {
                    cerr << "theColumn of qualifier = " << theColumn << endl;
                    cerr << "line = " << line << endl;
                    cerr << "tokens.size() of line = " << tokens.size() << endl;
                    cerr << "theColumn is >= tokens.size().  Exiting." << endl;
                    exit(1);
                }
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



/** \file FeaturesTableReader.hpp
 * Read #Features from a file with no additional qualifiers.
 */
template<typename LocType>
void readFeaturesTable(Features<LocType>& theFeatures, string filename,
                       long idPos = 0,
                       size_t contigColumn = 1,
                       long strandColumn = 2,
                       size_t startColumn = 3,
                       size_t endColumn = 4,
                       bool sorted = false,
                       string delimiter = "\t",
                       bool haveHeader = true,
                       bool trimWhiteSpaces = false,
                       bool collapseToStart = false,
                       bool collapseToEnd = false,
                       typename LocType::position_type::element_type startShift = 0,
                       typename LocType::position_type::element_type endShift = 0,
                       size_t resLength = 0,
                       string defaultStrand = "+",
                       long fileIndex = 0,
                       string senseStrand = "+" )
{
#ifdef DEBUG
    assert((collapseToStart && collapseToEnd) == false);
#endif
    StringQualifiers noQualifiers;
    readFeaturesTable(theFeatures, filename, noQualifiers,
                      idPos, contigColumn, strandColumn, startColumn, endColumn,
                      sorted, delimiter, haveHeader, trimWhiteSpaces,
                      collapseToStart, collapseToEnd, startShift, endShift,
                      resLength, defaultStrand, fileIndex, senseStrand);
}

#endif

