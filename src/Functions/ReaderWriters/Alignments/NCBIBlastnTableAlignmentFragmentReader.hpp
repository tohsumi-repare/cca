#ifndef MOLBIOLIB_NCBIBLASTNTABLEALIGNMENTFRAGMENTREADER_H
#define MOLBIOLIB_NCBIBLASTNTABLEALIGNMENTFRAGMENTREADER_H 1

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
#include "src/Objects/Interval.hpp"
#include "src/Objects/Location.hpp"
#include "src/Functions/ReaderWriters/Alignments/SequenceAlignmentFragmentReader.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/SystemUtilities/FileModificationTime.hpp"
#include "src/Functions/SystemUtilities/Cstdio.hpp"



/** \file NCBIBlastnTableAlignmentFragmentReader.hpp
 * Defines a NCBI blastn table (-m 8) reader.
 */


/* Below is not to be included in the doxygen documentation,
 * so only one * is used. */
/* Holds NCBI Blastn alignment information.
 * The names of these correspond to entries found in NCBI blast documentation
 * such as http://kb.iu.edu/data/awvn.html.  Not to be used by users of
 * MolBioLib.
 */
class NCBIBlastnAlignmentFragmentInfo
{
public:
    double percentIdentity, eValue, bitScore;
    size_t alignmentLength, numMismatches, numGaps;
};





/** The NCBI Blastn tab-separated file format with no header (-m 8 option) reader.
 * Derived from #SequenceAlignmentFragmentReader.  For usage,
 * see the documentation and comments of #SequenceAlignmentFragmentReader.
 * Exactly the same usage as
 * #BroadQltoutAlignmentFragmentReader, except use #NCBIBlastnTableAlignmentFragmentReader and
 * a blastn tabbed output file using option -m 8.  Subject sequence means the
 * matching sequence in the database.
 *
 * The method numErrors() (of type <code>unsigned long</code>) has been added
 * to this class, which includes both SNPs and indels.
 *
 * The alignmentScore() returns the eValue score.
 *
 * \todo Add functionality to theTargetIntervalIndex.
 */
class NCBIBlastnTableAlignmentFragmentReader : public SequenceAlignmentFragmentReader
{
public:
    // Do not need to override open() method
    // from SequenceAlignmentFragmentReader.

    NCBIBlastnTableAlignmentFragmentReader(bool theDoQueryIndex = true, bool theDoTargetIndex = true, bool theDoLocationIndex = true, bool theSequentialMode = false, float maxLoadFactor = 1.0)
    {
        doQueryIndex = theDoQueryIndex;
        doTargetIndex = theDoTargetIndex;
        doLocationIndex = theDoLocationIndex;
        sequentialMode = theSequentialMode;
        string targetLocExt = "";
        string queryExt = "", targetExt = "";
        limitReads = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        this->setIndexFileExtensions(queryExt, targetExt);
        setTargetLocationIndexFileExtension(targetLocExt);
        this->setMaxLoadFactor(maxLoadFactor);
    }

    NCBIBlastnTableAlignmentFragmentReader(string filename, string queryExt = "", string targetExt = "", string targetLocExt = "", bool theDoQueryIndex = true, bool theDoTargetIndex = true, bool theDoLocationIndex = true, bool theSequentialMode = false, float maxLoadFactor = 1.0)
    {
        doQueryIndex = theDoQueryIndex;
        doTargetIndex = theDoTargetIndex;
        doLocationIndex = theDoLocationIndex;
        sequentialMode = theSequentialMode;
        limitReads = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        this->setIndexFileExtensions(queryExt, targetExt);
        setTargetLocationIndexFileExtension(targetLocExt);
        this->setMaxLoadFactor(maxLoadFactor);
        open(filename);
    }


    size_t numErrors()
    {
        size_t numMismatches = convertFromString<size_t>((((this->currentFragment).getQueryLocation()).info()).getQualifier("numMismatches"));
        size_t numGaps = convertFromString<size_t>((((this->currentFragment).getQueryLocation()).info()).getQualifier("numGaps"));
        return numMismatches + numGaps;
    }


    double alignmentScore()
    {
        return convertFromString<size_t>((((this->currentFragment).getQueryLocation()).info()).getQualifier("eValue"));
    }


    bool readAlignmentFragment(size_t n)
    {
        // Ignore readerMode flag.  This is fast enough as is.

        if (n >= theFileTable.size())
        {
            return false;
        }

        currentFragment.clear();

        tempLine = theFileTable[n];
        splitString(tempLine, "\t", tempTokens);

        theOrigin = tempQueryLocation.origin();

        tempQueryLocation.contig() = tempTokens[0];
        tempSubjectLocation.contig() = tempTokens[1];
        tempInfoQuery.clear();
        tempInfoSubject.clear();
        tempInfoQuery.push_back("percentIdentity", tempTokens[2]);
        tempInfoSubject.push_back("percentIdentity", tempTokens[2]);
        tempInfoQuery.push_back("alignmentLength", tempTokens[3]);
        tempInfoSubject.push_back("alignmentLength", tempTokens[3]);
        tempInfoQuery.push_back("numMismatches", tempTokens[4]);
        tempInfoSubject.push_back("numMismatches", tempTokens[4]);
        tempInfoQuery.push_back("numGaps", tempTokens[5]);
        tempInfoSubject.push_back("numGaps", tempTokens[5]);
        start = convertFromString<long>(tempTokens[6]);
        end = convertFromString<long>(tempTokens[7]);
        if (start <= end)
        {
            tempSubjectLocation.strand() = "+";
        }
        else
        {
            tempSubjectLocation.strand() = "-";
            swap(start, end);
        }
        tempIntervalQuery.assign(start - 1 + theOrigin, end - 1 + theOrigin);
        start = convertFromString<long>(tempTokens[8]);
        end = convertFromString<long>(tempTokens[9]);
#ifdef DEBUG
        assert(start <= end);   // Currently can only handle RC on query.
#endif
        tempQueryLocation.strand() = "+";
        tempIntervalSubject.assign(start - 1 + theOrigin, end - 1 + theOrigin);
        tempInfoQuery.push_back("eValue", tempTokens[10]);
        tempInfoSubject.push_back("eValue", tempTokens[10]);
        tempInfoQuery.push_back("bitScore", tempTokens[11]);
        tempInfoSubject.push_back("bitScore", tempTokens[11]);

        if (tempInfoQuery.getQualifier("numGaps") == "0")
        {
            if (tempInfoQuery.getQualifier("numMismatches") == "0")
            {
                anElement.type() = NoChange;
            }
            else
            {
                anElement.type() = SNPs;
            }
        }
        else
        {
            if (tempInfoQuery.getQualifier("numMismatches") == "0")
            {
                anElement.type() = Indel;
            }
            else
            {
                anElement.type() = IndelSNPs;
            }
        }

        tempQueryLocation.position() = tempIntervalQuery;
        tempQueryLocation.info() = tempInfoQuery;
        currentFragment.globalQueryLocation = tempQueryLocation;
        anElement.queryLocation() = tempQueryLocation;

        tempSubjectLocation.position() = tempIntervalSubject;
        tempSubjectLocation.info() = tempInfoSubject;
        anElement.targetLocation() = tempSubjectLocation;
        currentFragment.globalTargetLocation = tempSubjectLocation;

        currentFragment.push_back(anElement);

        return true;
    }


    void setTargetLocationIndexFileExtension(string targetLocExt = "")
    {
        if (targetLocExt == "")
        {
            targetLocationIndexFileExtension = "MolBioLib.NCBIBlastnTableAlignmentFragmentReader.TargetLocation.index";
        }
        else
        {
            targetLocationIndexFileExtension = targetLocExt;
        }
    }


    bool paired(NCBIBlastnTableAlignmentFragmentReader& rhs)
    {
        // NCBI BLast is never paired.
        return false;
    }


    bool rcQueries()
    {
        return true;
    }


    bool hasSequences()
    {
        return false;
    }


    bool hasQualities()
    {
        return false;
    }


protected:
    void createTargetLocationIndices(string filename)
    {
        // TO BE WRITTEN
        // Check that filename + extension does not exist or else is deleted.
        // Override this function in derived class.
    }


    void openTargetLocationIndices(string filename)
    {
        // TO BE WRITTEN
        // Override this function in derived class.
    }


    ifstream::pos_type checkTargetLocationFile(string filename)
    {
        string tempFileName = filename + "." + targetLocationIndexFileExtension;
#ifdef DEBUG
        cerr << "Until NCBIBlastnTableAlignmentFragmentReader::createTargetLocationIndices is written, NCBIBlastnTableAlignmentFragmentReader::checkTargetLocationFile will always return 1.\n";
        return 1;
#endif
        return fileSize(tempFileName);
    }




private:
    string tempLine;
    vector<string> tempTokens;
    SequenceLocation tempQueryLocation, tempSubjectLocation;
    SequenceLocation::element_type theOrigin;
    long start, end;
    SequenceLocation::position_type tempIntervalQuery, tempIntervalSubject;
    SequenceLocation::info_type tempInfoQuery, tempInfoSubject;
    SequenceAlignmentFragment::alignment_element_type anElement;

};

#endif

