#ifndef MOLBIOLIB_BROADQLTOUTALIGNMENTFRAGMENTREADER_H
#define MOLBIOLIB_BROADQLTOUTALIGNMENTFRAGMENTREADER_H 1

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


#include "src/Objects/Interval.hpp"
#include "src/Objects/Location.hpp"
#include "src/Objects/Qualifiers.hpp"
#include "src/Objects/ReadOnlyStringFile.hpp"
#include "src/Functions/ReaderWriters/Alignments/SequenceAlignmentFragmentReader.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"


/** \file BroadQltoutAlignmentFragmentReader.hpp
 * Broad QLTOUT reader
 * \def DEFAULT_CONTIG_STRING
 * Default contig string name to put in front of query number.
 */
#define DEFAULT_CONTIG_STRING "chr"




/** \file BroadQltoutAlignmentFragmentReader.hpp
 * Broad QLTOUT reader
 * \enum BroadContigTranslation
 * \todo Replace <code>enum BroadContigTranslations</code> with
 * <code>enum class BroadContigTranslations</code> when strongly typed
 * enums work.
 */
enum BroadContigTranslation
{
    BroadContigNumber, BroadContigNumberDefaultChr, BroadContigStringQualifiers
};

/** The Broad qltout file format reader.
 * Derived from #SequenceAlignmentFragmentReader.  See #SequenceAlignmentFragmentReader
 * documentation and comments.  Note that we have added the method numErrors()
 * (of type <code>long</code>) to
 * this class that returns the number of SNPs and indels in the alignment.
 *
 * alignmentScore() returns the number of errors.
 *
 * Note that while the Broad's intervals are half-open, this reader will
 * convert them to closed intervals.  Also, the 0-based index is converted to
 * whatever index SequenceLocation's is.
 *
 * \todo Add functionality to theTargetIntervalIndex.
 */
class BroadQltoutAlignmentFragmentReader : public SequenceAlignmentFragmentReader
{
public:

    size_t numErrors()
    {
        StringQualifiers& theInfo = ((this->currentFragment).getTargetLocation()).info();
        return convertFromString<size_t>(theInfo.getQualifier("global_numErrors"));
    }

    double alignmentScore()
    {
        return static_cast<double>(numErrors());
    }

    void open(string filename)
    {
        // We need to override this function since we also need to know if this
        // qltout file is human readable or not and which line is the first QUERY
        // line.
        theFileName = filename;
        theFileTable.setSequentialMode(sequentialMode);
        theFileTable.open(filename);
        currentFileIndex = 0;
        // Find first QUERY line.
        while (currentFileIndex < theFileTable.size() &&
                theFileTable[currentFileIndex].substr(0, 5) != "QUERY")
        {
            ++currentFileIndex;
        }
        headerOffset = currentFileIndex;
        // Check if human-readable by checking if next line has QUERY in it.
        humanReadable = true;
        if (theFileTable.size() <= headerOffset + 1 ||
                theFileTable[currentFileIndex+1].substr(0, 5) == "QUERY")
        {
            humanReadable = false;
        }

        // Check if last line has done date in it.
        if (theFileTable.size() > 0)
        {
            if (theFileTable[theFileTable.size()-1].find(": done.") != string::npos)
            {
                lastLineDoneDate = true;
            }
            else
            {
                lastLineDoneDate = false;
            }
        }
        else
        {
            lastLineDoneDate = false;
        }

        currentFileIndex = -1;

        this->openIndices();
    }

    BroadQltoutAlignmentFragmentReader(BroadContigTranslation whichQueryTrans = BroadContigNumber, StringQualifiers& theQueryQualifiers = nullStringQualifiers, BroadContigTranslation whichTargetTrans = BroadContigNumber, StringQualifiers& theTargetQualifiers = nullStringQualifiers, bool theDoQueryIndex = true, bool theDoTargetIndex = true, bool theDoLocationIndex = true, bool theSequentialMode = false, float maxLoadFactor = 1.0)
    {
        ptrQueryTranslationTable = &theQueryQualifiers;
        ptrTargetTranslationTable = &theTargetQualifiers;
        doQueryIndex = theDoQueryIndex;
        doTargetIndex = theDoTargetIndex;
        doLocationIndex = theDoLocationIndex;
        sequentialMode = theSequentialMode;
        theQueryTranslation = whichQueryTrans;
        theTargetTranslation = whichTargetTrans;
        string queryExt = "", targetExt = "", indexExt ="", targetLocExt = "";
        this->setIndexFileExtensions(queryExt, targetExt, indexExt);
        setTargetLocationIndexFileExtension(targetLocExt);
        humanReadable = false;
        limitReads = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        this->setMaxLoadFactor(maxLoadFactor);
    }

    BroadQltoutAlignmentFragmentReader(string filename, BroadContigTranslation whichQueryTrans = BroadContigNumber, StringQualifiers& theQueryQualifiers = nullStringQualifiers, BroadContigTranslation whichTargetTrans = BroadContigNumber, StringQualifiers& theTargetQualifiers = nullStringQualifiers, string queryExt = "", string targetExt = "", string indexExt = "", string targetLocExt = "", bool theDoQueryIndex = true, bool theDoTargetIndex = true, bool theDoLocationIndex = true, bool theSequentialMode = false, float maxLoadFactor = 1.0)
    {
        ptrQueryTranslationTable = &theQueryQualifiers;
        ptrTargetTranslationTable = &theTargetQualifiers;
        doQueryIndex = theDoQueryIndex;
        doTargetIndex = theDoTargetIndex;
        doLocationIndex = theDoLocationIndex;
        sequentialMode = theSequentialMode;
        theQueryTranslation = whichQueryTrans;
        theTargetTranslation = whichTargetTrans;
        this->setIndexFileExtensions(queryExt, targetExt, indexExt);
        setTargetLocationIndexFileExtension(targetLocExt);
        this->setMaxLoadFactor(maxLoadFactor);
        limitReads = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        open(filename);
    }

    // Below is probably what will be used most often.
    BroadQltoutAlignmentFragmentReader(string filename, StringQualifiers& theTargetQualifiers, bool theDoQueryIndex = true, bool theDoTargetIndex = true, bool theDoLocationIndex = true, bool theSequentialMode = false, float maxLoadFactor = 1.0)
    {
        ptrQueryTranslationTable = &nullStringQualifiers;
        ptrTargetTranslationTable = &theTargetQualifiers;
        doQueryIndex = theDoQueryIndex;
        doTargetIndex = theDoTargetIndex;
        doLocationIndex = theDoLocationIndex;
        sequentialMode = theSequentialMode;
        theQueryTranslation = BroadContigNumber;
        theTargetTranslation = BroadContigStringQualifiers;
        string queryExt = "", targetExt = "", indexExt = "", targetLocExt = "";
        this->setIndexFileExtensions(queryExt, targetExt, indexExt);
        setTargetLocationIndexFileExtension(targetLocExt);
        this->setMaxLoadFactor(maxLoadFactor);
        limitReads = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        open(filename);
    }

    size_t size()
    {
        size_t numEntries = theFileTable.size() - headerOffset;
        if (lastLineDoneDate && (numEntries > 0))
        {
            --numEntries;
        }
        if (humanReadable)
        {
            numEntries /= 2;
        }
        return numEntries;
    }


    bool readAlignmentFragment(size_t n)
    {

        currentFragment.clear();

        tempIndex = headerOffset;
        if (humanReadable)
        {
            tempIndex += 2*n;
        }
        else
        {
            tempIndex += n;
        }

        // -1 because the last line is usually the done date.  :-(
        if (lastLineDoneDate)
        {
            if (tempIndex >= theFileTable.size()-1)
            {
                return false;
            }
        }
        else
        {
            if (tempIndex >= theFileTable.size())
            {
                return false;
            }
        }


        tempLine = theFileTable[tempIndex];
        splitString(tempLine, "\t", tempTokens);

#ifdef DEBUG
        if (tempTokens[0] != "QUERY")
        {
            cerr << "Error in BroadQltoutAlignmentFragmentReader!  The alignment file must contain only a header, then just QUERY lines (with associated human lines, if human-readable) only.  Another header cannot be placed in the alignment file.  Exiting.\n";
            assert(true == false);
        }
#endif

        // Now that we have the tokens, can fill currentFragment.
        // Recall stop is +1 beyond the final position.
        q = tempTokens[1];
        // currentQ is used so that we can easily tell if this alignment is
        // paired with another.
        currentQ = q;
        // Broad index starts at 0 so adjust by origin
        // E.g. if the origin is 1, then add 1 to all Broad coordinates.
        long seqLocOrigin = tempQueryLocation.origin();
        q_start = convertFromString<long>(tempTokens[2]) + seqLocOrigin;
        q_stop = convertFromString<long>(tempTokens[3]) + seqLocOrigin;
        q_length = convertFromString<long>(tempTokens[4]);
#ifdef PROG_DEBUG
        // We cannot yet handle the case where not all of it aligns.
        // Else must load in reads and ref and check for
        // beginning and ending indels.
        assert(q_length == q_stop - q_start);
#endif
        q_rc = tempTokens[5];
        t = tempTokens[6];
        // Broad index starts at 0.
        t_start = convertFromString<long>(tempTokens[7]) + seqLocOrigin;
        t_stop = convertFromString<long>(tempTokens[8]) + seqLocOrigin;
        t_length = tempTokens[9];
        n_blocks = convertFromString<size_t>(tempTokens[10]);

        tempInterval.assign(q_start, q_stop - 1);  // Broad uses half-open
        // intervals.
        tempQueryLocation.strand() = "+";
        size_t theQueryNum;
        switch (theQueryTranslation)
        {
            case BroadContigNumber :
                tempQueryLocation.contig() = q;
                break;
            case BroadContigNumberDefaultChr :
                theQueryNum = convertFromString<size_t>(q) + 1;
                tempQueryLocation.contig() = DEFAULT_CONTIG_STRING +
                                             convertToString<size_t>(theQueryNum);
                break;
            case BroadContigStringQualifiers :
                tempQueryLocation.contig() = ptrQueryTranslationTable->getQualifier(q);
                break;
        }
        tempQueryLocation.position() = tempInterval;
        currentFragment.globalQueryLocation = tempQueryLocation;
        currentFragment.globalQueryLocation.info().push_back("query_contig_length",
                convertToString<long>(q_length));
        // globalQueryLocation.info() will be set after blocks (tot errors)
        tempInterval.assign(t_start, t_stop - 1);
        size_t theTargetNum;
        switch (theTargetTranslation)
        {
            case BroadContigNumber :
                tempTargetLocation.contig() = t;
                break;
            case BroadContigNumberDefaultChr :
                theTargetNum = convertFromString<size_t>(t) + 1;
                tempTargetLocation.contig() = DEFAULT_CONTIG_STRING +
                                              convertToString<size_t>(theTargetNum);
                break;
            case BroadContigStringQualifiers :
                tempTargetLocation.contig() = ptrTargetTranslationTable->getQualifier(t);
                break;
        }
        if (q_rc == "0")
        {
            tempTargetLocation.strand() = "+";
        }
        else
        {
            tempTargetLocation.strand() = "-";
        }
        tempTargetLocation.position() = tempInterval;
        currentFragment.globalTargetLocation = tempTargetLocation;
        // globalTargetLocation.info()'s total error will be set after blocks
        // Need below to reconstruct orignal qltout when writing.
        currentFragment.globalTargetLocation.info().push_back("target_contig_length",
                t_length);

        long numErrors = 0;


        if (readerMode == NoPoly)
        {
            n_blocks = 1;
            g = 0;
            e = 0;
            b = q_length;
            tempQueryLocation.info().clear();
            tempTargetLocation.info().clear();
            anElement.type() = NoChange;
            tempTargetLocation.info().push_back("element_gapLength", "0");
            q_end = q_start + b - 1;
            t_end = t_start + b - 1;
            tempInterval.assign(q_start, q_end);
            tempQueryLocation.position() = tempInterval;
            anElement.queryLocation() = tempQueryLocation;
            tempInterval.assign(t_start, t_end);
            tempTargetLocation.position() = tempInterval;
            anElement.targetLocation() = tempTargetLocation;
            currentFragment.push_back(anElement);
            // No need to advance, since only one block pushed.

        }
        else
        {
            // readerMode == Normal
            for (size_t i = 0; i < n_blocks; ++i)
            {
                g = convertFromString<long>(tempTokens[11 + 3*i]);
                b = convertFromString<size_t>(tempTokens[11 + 3*i + 1]);
                e = convertFromString<size_t>(tempTokens[11 + 3*i + 2]);

                tempQueryLocation.info().clear();
                tempTargetLocation.info().clear();

                // Start at q_start and t_start


                if (e == 0)
                {
                    // Add alignment row
                    anElement.type() = NoChange;
                }
                else
                {
                    // Add SNPs row
                    anElement.type() = MatchOrMismatch;
                    // SNP errors.  Don't know where they occur.
                    numErrors += e;
                }
                tempTargetLocation.info().push_back("element_numErrors", convertToString<size_t>(e));

                if (g > 0)
                {
                    anElement.type() = Deletion;
                    q_end = q_start;
                    t_end = t_start + g - 1;

                }
                else if (g < 0)
                {
                    anElement.type() = Insertion;
                    q_end = q_start - g - 1;
                    t_end = t_start;
                }
                if (g != 0)
                {
                    tempInterval.assign(q_start, q_end);
                    tempQueryLocation.position() = tempInterval;
                    long absG = (g > 0) ? g : -g;
                    tempTargetLocation.info().push_back("element_gapLength", convertToString<long>(absG));
                    if ((anElement.type() != Deletion) ||
                            (anElement.type() != Insertion))
                    {
                        numErrors += static_cast<size_t>(absG);
                    }
                    anElement.queryLocation() = tempQueryLocation;

                    tempInterval.assign(t_start, t_end);
                    tempTargetLocation.position() = tempInterval;
                    anElement.targetLocation() = tempTargetLocation;

                    if (g > 0)
                    {
                        t_start += absG;
                    }
                    else if (g < 0)
                    {
                        q_start += absG;
                    }

                    currentFragment.push_back(anElement);
                    // Advance positions to after last interval end.
                    if (g < 0)
                    {
                        q_start = q_end + 1;
                    }
                    if (g > 0)
                    {
                        t_start = t_end + 1;
                    }

                    anElement.type() = NoChange;
                }

                if (b > 0)
                {
                    // tempTargetLocation.info().push_back("element_gapLength", "0");
                    q_end = q_start + b - 1;
                    t_end = t_start + b - 1;
                    tempInterval.assign(q_start, q_end);
                    tempQueryLocation.position() = tempInterval;
                    anElement.queryLocation() = tempQueryLocation;
                    tempInterval.assign(t_start, t_end);
                    tempTargetLocation.position() = tempInterval;
                    anElement.targetLocation() = tempTargetLocation;
                    currentFragment.push_back(anElement);
                    // Advance positions to after last interval end.
                    q_start = q_end + 1;
                    t_start = t_end + 1;
                }
            }
        }

        currentFragment.globalTargetLocation.info().push_back("global_numErrors", convertToString<long>(numErrors));

        return true;
    }



    void setIndexFileExtensions(string queryExt = "", string targetExt = "", string extensionExt = "")
    {
        string queryContigTranslationString = convertBroadContigTranslationToString(theQueryTranslation);
        string targetContigTranslationString = convertBroadContigTranslationToString(theTargetTranslation);
        if (queryExt == "")
        {
            queryContigIndexFileExtension = "MolBioLib.BroadQltoutAlignmentFragmentReader." + queryContigTranslationString + ".QueryContig.index";
        }
        else
        {
            queryContigIndexFileExtension = queryExt;
        }
        if (targetExt == "")
        {
            targetContigIndexFileExtension = "MolBioLib.BroadQltoutAlignmentFragmentReader." + targetContigTranslationString + ".TargetContig.index";
        }
        else
        {
            targetContigIndexFileExtension = targetExt;
        }
        if (extensionExt == "")
            // Below is so that the final extension is always .index
        {
            indexFileExtension = "IndexOf.index";
        }
        else
        {
            indexFileExtension = extensionExt;
        }
    }



    void setTargetLocationIndexFileExtension(string targetLocExt = "")
    {
        if (targetLocExt == "")
        {
            targetLocationIndexFileExtension = "MolBioLib.BroadQltoutAlignmentFragmentReader.TargetLocation.index";
        }
        else
        {
            targetLocationIndexFileExtension = targetLocExt;
        }
    }


    bool paired(BroadQltoutAlignmentFragmentReader& rhs)
    {
        if (currentQ == rhs.currentQ)
        {
            return true;
        }
        else
        {
            return false;
        }
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

    // string getCurrentSequence() and string getCurrentQuality() are not
    // available.  Also, Fasta& getSequences() and Qual& getQualities() are not
    // available either.  These are because the qltout format does not hold
    // the original quality or sequence.

    ~BroadQltoutAlignmentFragmentReader()
    {
        tempQueryLocation.info().clear();
        tempQueryLocation.info().clear();
        tempTargetLocation.strand().clear();
        tempTargetLocation.strand().clear();
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
        cerr << "Until BroadQltoutAlignmentFragmentReader::createTargetLocationIndices is written, BroadQltoutAlignmentFragmentReader::checkTargetLocationFile will always return 1.\n";
        return 1;
#endif
        return fileSize(tempFileName);
    }


private:
    BroadContigTranslation theQueryTranslation, theTargetTranslation;
    StringQualifiers *ptrQueryTranslationTable, *ptrTargetTranslationTable;
    bool humanReadable;
    bool lastLineDoneDate;
    size_t headerOffset, tempIndex;
    string tempLine;
    vector<string> tempTokens;
    string q, t, q_rc;
    long q_start, q_stop, q_length, t_start, t_stop;
    string t_length;
    size_t n_blocks;
    long g, b, e;
    long q_end, t_end;
    SequenceLocation tempQueryLocation, tempTargetLocation;
    SequenceLocation::position_type tempInterval;
    string tempContig;
    SequenceAlignmentFragment::alignment_element_type anElement;
    string currentQ;


    string convertBroadContigTranslationToString(BroadContigTranslation x)
    {
        switch(x)
        {
            case BroadContigNumber :
                return "BCN";
                break;
            case BroadContigNumberDefaultChr :
                return "BCNDC";
                break;
            case BroadContigStringQualifiers :
                return "BCSQ";
                break;
            default:
                return "NA" ;
                break;
        }
    }


};



/** \file BroadQltoutAlignmentFragmentReader.hpp
 * \fn void createBroadMapStringQualifiersFromLabelsFile(string filename, StringQualifiers& mapChr)
 * Given a labels file (just the names without the &gt; of a
 * reference fasta file),
 * this generates the mapping needed for the alignment reader.
 */
void createBroadMapStringQualifiersFromLabelsFile(string filename,
        StringQualifiers& mapChr)
{
    string line;
    ReadOnlyStringFile labelsTableFile(filename, "createBroadMapStringQualifiersFromLabelsFile.index", false);
    size_t labelsTableSize = labelsTableFile.size();
    string currI;
    for (size_t i = 0; i < labelsTableSize; ++i)
    {
        currI = convertToString<size_t>(i);
        line = labelsTableFile[i];
        mapChr.push_back(currI, line);
    }
    labelsTableFile.close();
}



/** \file BroadQltoutAlignmentFragmentReader.hpp
 * \fn void createBroadMapStringQualifiersFromHeaderSizes(string filename, StringQualifiers& mapChr)
 * Given a headerSizes file (header [tab] size [newline] of a reference fasta
 * or fastq file), this generates the mapping needed for the alignment reader.
 */
void createBroadMapStringQualifiersFromHeaderSizes(string filename,
        StringQualifiers& mapChr)
{
    string line;
    ReadOnlyStringFile headerSizesTableFile(filename, "createBroadMapStringQualifiersFromHeaderSizes.index", false);
    size_t headerSizesTableSize = headerSizesTableFile.size();
    string currI;
    vector<string> tokens;
    for (size_t i = 0; i < headerSizesTableSize; ++i)
    {
        currI = convertToString<size_t>(i);
        line = headerSizesTableFile[i];
        splitString(line, "\t", tokens);
        mapChr.push_back(currI, tokens[0]);
    }
    headerSizesTableFile.close();
}




/** \file BroadQltoutAlignmentFragmentReader.hpp
 * \fn void createBroadMapStringQualifiersFromFastaFile(string filename, StringQualifiers& mapChr)
 * Given a FASTA or FASTQ file, this generates the mapping needed for the
 * alignment reader.
 */
void createBroadMapStringQualifiersFromFastaFile(string filename,
        StringQualifiers& mapChr)
{
    Fasta theReads(filename);
    size_t numReads = theReads.size();

    string currI;
    for (size_t i = 0; i < numReads; ++i)
    {
        currI = convertToString<size_t>(i);
        mapChr.push_back(currI, theReads.getHeader(i));
    }

    // theReads will be destroyed once it gets out of this scope.
}



/** \file BroadQltoutAlignmentFragmentReader.hpp
 * \fn void createBroadMapStringQualifiersFromFasta(Fasta& filename, StringQualifiers& mapChr)
 * Given a FASTA data that has been filled, this generates the mapping needed
 * for the alignment reader.
 */
void createBroadMapStringQualifiersFromFasta(Fasta& theReads,
        StringQualifiers& mapChr)
{
    size_t numReads = theReads.size();

#ifdef DEBUG
    if (numReads == 0)
    {
        cerr << "Warning from createBroadMapStringQualifiersFromFasta!  The Fasta  data structure passed has no reads in them.  Continuing execution." << endl;
    }
#endif

    string currI = "";
    for (size_t i = 0; i < numReads; ++i)
    {
        currI = convertToString<size_t>(i);
        mapChr.push_back(currI, theReads.getHeader(i));
    }
}

#endif



