#ifndef MOLBIOLIB_FEATURESTABLEENSEMBLSTRUCTURESREADER_H
#define MOLBIOLIB_FEATURESTABLEENSEMBLSTRUCTURESREADER_H 1

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
#include "src/Functions/ReaderWriters/Features/FeaturesTableEnsemblQualifiersReader.hpp"


/** \file FeaturesTableEnsemblStructuresReader.hpp
 * Read #Features from a mart_export.txt file from Ensembl's Biomart.
 * \fn void readFeaturesTableEnsemblStructures(SequenceFeatures& theFeatures, string featureFile = "mart_export.txt", bool sorted = false, bool collapseToStart = false, bool collapseToEnd = false, SequenceFeatures::element_type startShift = 0, SequenceFeatures::element_type endShift = 0, bool suppressExons = false, bool suppressIntrons = false, size_t resLength = 0)
 * Select Structures.  Then, you must
 * select <em>all</em> of the Feature attributes in order - first the left column,
 * then the right.  Then, select all of the Exon information - again first the
 * left column, then the right.
 * Example usage is:
 * \code
 * SequenceFeatures newFeatures;
 * readFeaturesTableEnsemblStructures(newFeatures,
 *                                    "mart_export.hg18.txt", false,
 *                                    false, false, 0, 0, false, false, 0);
 * \endcode
 * The last seven parameters are optional.  By default, the feature file is
 * "mart_export.txt".  mart_export.txt is not sorted, so the default is
 * correct.  However, if set to <code>true</code>, then #Features can search for
 * features hit at a particular #Location much faster.  The -1 is the reserve
 * length of the string used to read in one line of mart_export.txt.  -1 means
 * not to reserve any space and let the string class handle it.  Perhaps 2000 is
 * not a bad value here.  The 0 and 0 are the start and stop shifts of the
 * global feature coordinates only, like in #readFeaturesTable.  Exon
 * coordinates are not shifted.  The second
 * and third <code>false</code> specifies the position is not collapsed to the
 * start and end coordinates respectively.  Exon coordinates are not collapsed.
 * The fourth <code>false</code> says we
 * are not suppressing the creation of the exons rows in the #Features table.
 * The last <code>false</code> specifies we are not suppressing the creation of
 * intron rows in the #Features table.  Note that the
 * featureFile <em>must</em> be
 * processed by preprocessEnsemblStructuresFile.cpp first.
 *
 * The idea here is that we create one feature whose coordinates is the entire
 * feature and then for each row (i.e. each exon), we create another feature whose
 * coordinates is just that exon's.  Thus, we can see who hit the exons and
 * who hit the features (but not the exons).  The qualifier to see if a feature is a
 * complete feature is labeled "Feature", which is 1 if so and 0 otherwise.
 */
void readFeaturesTableEnsemblStructures(SequenceFeatures& theFeatures,
                                        string featureFile = "mart_export.txt",
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

    // Ensembl has a header.
    ReadOnlyStringFile ensemblTableFile(featureFile, true, true, true, "readFeaturesTableEnsemblStructures.index", resLength);
    for (size_t i = 0; !ensemblTableFile.fail(); ++i)
    {
        line = ensemblTableFile[i];
        splitString(line, "\t", tokens);
#ifdef DEBUG
        // Below only true if preprocessed.
        assert(tokens.size() == 29);
#endif

        SequenceLocation newLoc;
        SequenceLocation::element_type theOrigin = newLoc.origin();
        newLoc.contig() = tokens[3];
        if (tokens[8] == "1")
        {
            newLoc.strand() = "+";
        }
        else if (tokens[8] == "-1")
        {
            newLoc.strand() = "-";
        }
        SequenceLocation::position_type newI;
        // These are the actual starting coordinates of the object in row i.
        SequenceFeatures::element_type startPos = convertFromString<SequenceFeatures::element_type>(tokens[27]) - 1 + theOrigin;
        SequenceFeatures::element_type endPos = convertFromString<SequenceFeatures::element_type>(tokens[28]) - 1 + theOrigin;
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
        string currentFeatureID = tokens[0];
        newQualifiers.push_back("Ensembl Feature ID", currentFeatureID);
        newQualifiers.push_back("Ensembl Transcript ID", tokens[1]);
        newQualifiers.push_back("Ensembl Protein ID", tokens[2]);
        newQualifiers.push_back("Transcript Start (bp)", tokens[6]);
        newQualifiers.push_back("Transcript End (bp)", tokens[7]);
        newQualifiers.push_back("Associated Feature Name", tokens[9]);
        newQualifiers.push_back("Associated Feature DB", tokens[10]);
        newQualifiers.push_back("5' UTR Start", tokens[11]);
        newQualifiers.push_back("5' UTR End", tokens[12]);
        newQualifiers.push_back("3' UTR Start", tokens[13]);
        newQualifiers.push_back("3' UTR End", tokens[14]);
        newQualifiers.push_back("CDS Start", tokens[15]);
        newQualifiers.push_back("CDS End", tokens[16]);
        newQualifiers.push_back("CDS Length", tokens[17]);
        newQualifiers.push_back("Transcript count", tokens[18]);
        newQualifiers.push_back("Description", tokens[19]);
        newQualifiers.push_back("Feature Biotype", tokens[20]);
        newQualifiers.push_back("Ensembl Exon ID", tokens[21]);
        newQualifiers.push_back("Exon Chr Start (bp)", tokens[22]);
        newQualifiers.push_back("Exon Chr End (bp)", tokens[23]);
        newQualifiers.push_back("Exon Rank in Transcript", tokens[24]);
        newQualifiers.push_back("phase", tokens[25]);
        newQualifiers.push_back("Gene Section", tokens[26]);
        newQualifiers.push_back("Starting Coordinates of this Object", tokens[27]);
        newQualifiers.push_back("Ending Coordinates of this Object", tokens[28]);
        string currentIndex = convertToString<size_t>(i);
        if ((!suppressExons || (tokens[26].find("Exon") != string::npos)) &&
                (!suppressIntrons || (tokens[26].find("Intron") != string::npos)))
        {
            theFeatures.addFeature(currentIndex, newLoc, newQualifiers);
        }
    }
    ensemblTableFile.close();
}

#endif

