#ifndef MOLBIOLIB_HELICOSSMSALIGNMENTFRAGMENTREADER_H
#define MOLBIOLIB_HELICOSSMSALIGNMENTFRAGMENTREADER_H 1

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
#include "src/Objects/Qualifiers.hpp"
#include "src/Objects/Fasta.hpp"
#include "src/Objects/Qual.hpp"
#include "src/Objects/ReadOnlyStringFile.hpp"
#include "src/Functions/ReaderWriters/Alignments/SequenceAlignmentFragmentReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/BroadQltoutAlignmentFragmentReader.hpp"
#include "src/Functions/Transformers/StringOperations/ReplaceString.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"
#include "src/Functions/SystemUtilities/FileModificationTime.hpp"



/** \file HelicosSmsAlignmentFragmentReader.hpp
 * The Helicos sms format reader.
 * Derived from #SequenceAlignmentFragmentReader.  For usage, see
 * #SequenceAlignmentFragmentReader's
 * documentation/comments.  Note that we have added the method numErrors()
 * (of type <code>long</code>) to
 * this class that returns the number of SNPs and indels in the alignment.
 *
 * <code>alignmentScore()</code> returns MAPQ.
 *
 * <code>resetToBegin(true, false)</code> resets the sequential reading back
 * to the beginning of the file.  The true says to read the first part
 * of a paired read (if paired).  The false says to delete all the reads and
 * query name.  Those are the default values.
 *
 * \todo Add functionality to theTargetIntervalIndex.  Also work on the
 * paired() method.
 *
 * Specific to this reader, there is the function
 * <code>setRegernateHeaderMapsOnInput(bool newValue)</code>.  The value is
 * by default true, but suppose one is going to read in all the alignments
 * before accessing any of the sequences, then one can get the reads
 * and qualities and call regenerateHeaderMaps() on each afterwards.  In such a
 * case, one may want to pass false to the above function, so that the
 * readAlignmentFragment runs faster.
 */
class HelicosSmsAlignmentFragmentReader : public SequenceAlignmentFragmentReader
{
public:

    bool aligned()
    {
        // All alignment lines are aligned reads.
        return true;
    }


    vector<string> getHeaderSection()
    {
        return theHeaderSection;
    }

    bool rcQueries()
    {
        return false;
    }

    bool hasSequences()
    {
        return true;
    }

    bool hasQualities()
    {
        return false;
    }

    Fasta& getSequences()
    {
        return (*ptrSequences);
    }

    void setSequences(Fasta& theSeq)
    {
        ptrSequences = &theSeq;
    }

    void clearSequences()
    {
        ptrSequences->clear();
    }

    string getCurrentSequence()
    {
        return tagSeq;
    }
    string getCurrentReference()
    {
        return refSeq;
    }



    void open(string filename, string refHeaderSizesFile = "")
    {
        // We need to override this function since we also need to know if this
        // file has any header section.
        refHeaderSizes.clear();
        if (refHeaderSizesFile != "")
        {
            hasRefHeaderSizes = true;
            createBroadMapStringQualifiersFromHeaderSizes(refHeaderSizesFile, mapChr);
            refHeaderSizes.read(refHeaderSizesFile);

        }
        else
        {
            hasRefHeaderSizes = false;
        }
        theFileName = filename;
        theFileTable.setSequentialMode(sequentialMode);
        theFileTable.open(filename);
        currentFileIndex = 0;
        ptrSequences->clear();
        theHeaderSection.clear();
        // Find first line that does not start with # or Reference_ID.
        while ((currentFileIndex < theFileTable.size()) &&
                ((theFileTable[currentFileIndex].substr(0, 1) == "#") ||
                 (theFileTable[currentFileIndex].substr(0,12) == "Reference_ID")))
        {
            theHeaderSection.push_back(theFileTable[currentFileIndex]);
            ++currentFileIndex;
        }
        headerOffset = currentFileIndex;
        currentFileIndex = numeric_limits<size_t>::max();
        this->openIndices();
    }


    void resetToBegin(bool keepReads = false)
    {
        currentFileIndex = numeric_limits<size_t>::max();
        // Assumed the header is already read in.
        if (!keepReads)
        {
            ptrSequences->clear();
        }
    }


    // firstOrSecond = true -> do first read.
    HelicosSmsAlignmentFragmentReader(bool theKeepReads = true, bool theDoQueryIndex = true, bool theDoTargetIndex = true, bool theDoLocationIndex = true, bool theSequentialMode = false, float maxLoadFactor = 1.0)
    {
        keepReads = theKeepReads;
        doQueryIndex = theDoQueryIndex;
        doTargetIndex = theDoTargetIndex;
        doLocationIndex = theDoLocationIndex;
        sequentialMode = theSequentialMode;
        hasRefHeaderSizes = false;
        string queryExt = "", targetExt = "", indexExt = "", targetLocExt = "";
        limitReads = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        this->setIndexFileExtensions(queryExt, targetExt, indexExt);
        setTargetLocationIndexFileExtension(targetLocExt);
        this->setMaxLoadFactor(maxLoadFactor);
        ptrSequences = new Fasta;
        ptrSequences->setMaxLoadFactor(maxLoadFactor);
        regenerateHeaderMapsOnInsert = true;
        readSequences = true;
#ifdef DEBUG
        if (keepReads)
        {
            cerr << "Warning from HelicosSmsAlignmentFragmentReader!  Keeping the reads will make reading the SAM file much slower since it has to regenerate the header mappings on every read.  If not doing setReadSequences(false), consider batching and doing setRegenerateHeaderMapsOnInsert(false).  Continuing execution." << endl;
        }
#endif
    }

    HelicosSmsAlignmentFragmentReader(string filename, string refHeaderSizesFile = "", string queryExt = "", string targetExt = "", string indexExt = "", string targetLocExt = "", bool theKeepReads = true, bool theDoQueryIndex = true, bool theDoTargetIndex = true, bool theDoLocationIndex = true, bool theSequentialMode = false, float maxLoadFactor = 1.0)
    {
        keepReads = theKeepReads;
        doQueryIndex = theDoQueryIndex;
        doTargetIndex = theDoTargetIndex;
        doLocationIndex = theDoLocationIndex;
        sequentialMode = theSequentialMode;
        limitReads = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        this->setIndexFileExtensions(queryExt, targetExt, indexExt);
        setTargetLocationIndexFileExtension(targetLocExt);
        this->setMaxLoadFactor(maxLoadFactor);
        ptrSequences = new Fasta;
        ptrSequences->setMaxLoadFactor(maxLoadFactor);
        regenerateHeaderMapsOnInsert = true;
        readSequences = true;
#ifdef DEBUG
        if (keepReads)
        {
            cerr << "Warning from HelicosSmsAlignmentFragmentReader!  Keeping the reads will make reading the SAM file much slower since it has to regenerate the header mappings on every read.  If not doing setReadSequences(false), consider batching and doing setRegenerateHeaderMapsOnInsert(false).  Continuing execution." << endl;
        }
#endif
        open(filename, refHeaderSizesFile);
    }


    size_t size()
    {
        size_t numEntries = theFileTable.size() - headerOffset;
        return numEntries;
    }


    void regenerateHeaderMaps()
    {
        ptrSequences->regenerateHeaderMaps();
    }

    /* \fn size_t numErrors()
     * The problem is that with the CIGAR format, one cannot tell how many
     * SNPs are in [0-9]+M, since M means either match or mismatch.
     */
    size_t numErrors()
    {
        return numErrs;
    }

    bool readAlignmentFragment(size_t n)
    {

        currentFragment.clear();
        currentNumErrors = 0;

        tempIndex = headerOffset + n;

        if (tempIndex >= theFileTable.size())
        {
            return false;
        }

        line = theFileTable[tempIndex];
        splitString(line, "\t", tokens);

        isPaired = false;

        referenceID = tokens[0];
        tagID = tokens[1];
        referenceStart = convertFromString<long>(tokens[2]);
        referenceEnd = convertFromString<long>(tokens[3]);
        tagStart = convertFromString<long>(tokens[4]);
        tagEnd = convertFromString<long>(tokens[5]);
        nScore = convertFromString<double>(tokens[6]);
        numMatch = convertFromString<size_t>(tokens[7]);
        numDel = convertFromString<size_t>(tokens[8]);
        numIns = convertFromString<size_t>(tokens[9]);
        numSub = convertFromString<size_t>(tokens[10]);
        strand = tokens[11];
        tagSeq = tokens[12];
        refSeq = tokens[13];

        if (readerMode == NoPoly)
        {
            string temp = removeCharFromString(tagSeq, '-');
            tagSeq = temp;
            temp = removeCharFromString(refSeq, '-');
            refSeq = temp;
            size_t minLen = min(tagSeq.size(), refSeq.size());
            tagStart = 1;
            tokens[4] = "1";
            tagEnd = minLen;
            tokens[5] = convertToString<long>(tagEnd);
            referenceEnd = referenceStart + static_cast<long>(minLen) - 1;
            tokens[3] = convertToString<long>(referenceEnd);
            refSeq = refSeq.substr(0, minLen);
            tagSeq = refSeq;
            numMatch = minLen;
            tokens[7] = convertToString<size_t>(numMatch);
            numDel = 0;
            tokens[8] = "0";
            numIns = 0;
            tokens[9] = "0";
            numSub = 0;
            tokens[10] = "0";
        }

        tagSeqWithoutDashes = removeCharFromString(tagSeq, '-');

        if (readSequences && keepReads)
        {
            if (tagIDs.find(tagID) == tagIDs.end())
            {
                if (tagSeq != "")
                {
                    tagIDs.insert(tagID);
                    // Remove dashes, since they are now allowed in Fasta format.
                    ptrSequences->push_back(tagID, tagSeqWithoutDashes);
                    if (regenerateHeaderMapsOnInsert)
                    {
                        ptrSequences->regenerateHeaderMaps();
                    }
                }
                else
                {
                    string blankString = "";
                    Sequence tempSeq = blankString;
                    tagIDs.insert(tagID);
                    ptrSequences->push_back(tagID, tempSeq);
                    if (regenerateHeaderMapsOnInsert)
                    {
                        ptrSequences->regenerateHeaderMaps();
                    }
                }
            }
        }

        long seqLocOrigin = tempQueryLocation.origin();
        // -1 below since POS starts from 1.
        // From the examples, it seems that Helicos uses 1-based coordinates.
        tagStart -= 1;
        tagEnd -= 1;
        referenceStart += seqLocOrigin - 1;
        referenceEnd += seqLocOrigin - 1;

        tempGlobalQueryLocation.info().clear();
        tempGlobalQueryLocation.contig() = tagID;
        tempGlobalQueryLocation.strand() = "+";
        tempInterval.assign(tagStart, tagEnd);
        tempGlobalQueryLocation.position() = tempInterval;

        tempGlobalTargetLocation.info().clear();
        if (hasRefHeaderSizes)
        {
            // Subtract 1, mapChr starts from 0
            referenceID = convertToString<size_t>(convertFromString<size_t>(referenceID) - 1);
            referenceID = mapChr.getQualifier(referenceID);
        }
        tempGlobalTargetLocation.contig() = referenceID;
        tempGlobalTargetLocation.strand() = strand;
        tempInterval.assign(referenceStart, referenceEnd);
        tempGlobalTargetLocation.position() = tempInterval;
        tempGlobalTargetLocation.info().push_back("#_match", tokens[7]);
        tempGlobalTargetLocation.info().push_back("#_del", tokens[8]);
        tempGlobalTargetLocation.info().push_back("#_ins", tokens[9]);
        tempGlobalTargetLocation.info().push_back("#_sub", tokens[10]);
        tempGlobalTargetLocation.info().push_back("NScore", tokens[6]);
        numErrs = numDel + numIns + numSub;
        tempGlobalTargetLocation.info().push_back("global_numErrors", convertToString<long>(numErrs));

        if (readSequences)
        {
            tempGlobalQueryLocation.info().push_back("TagSeq", tagSeq);
            tempGlobalTargetLocation.info().push_back("RefSeq", refSeq);
        }

        tempGlobalQueryLocation.info().push_back("query_contig_length", convertToString<size_t>(tagSeqWithoutDashes.size()));
        if (hasRefHeaderSizes)
        {
            tempGlobalTargetLocation.info().push_back("target_contig_length", convertToString<size_t>(refHeaderSizes.getSize(referenceID)));
        }

        tempQueryLocation.contig() = tempGlobalQueryLocation.contig();
        tempQueryLocation.strand() = tempGlobalQueryLocation.strand();
        tempTargetLocation.contig() = tempGlobalTargetLocation.contig();
        tempTargetLocation.strand() = tempGlobalTargetLocation.strand();

        if (readerMode == NoPoly)
        {
            tempQueryLocation.info().clear();
            tempTargetLocation.info().clear();
            anElement.type() = NoChange;
            tempInterval.assign(tagStart, tagEnd);
            tempQueryLocation.position() = tempInterval;
            anElement.queryLocation() = tempGlobalQueryLocation;
            tempInterval.assign(referenceStart, referenceEnd);
            tempTargetLocation = tempGlobalTargetLocation;
            tempTargetLocation.position() = tempInterval;
            tempTargetLocation.info().push_back("element_gapLength", "0");
            anElement.targetLocation() = tempTargetLocation;
            currentFragment.push_back(anElement);
        }
        else
        {
            // First, variables that point to the location in tagSeq and refSeq:
            size_t currentQueryStart = static_cast<size_t>(tagStart);
            size_t currentQueryPos = currentQueryStart;
            // currentTargetStart never used.
            // size_t currentTargetStart = 0;
            size_t currentTargetPos = 0;
            // Next, variables that point to the location in tagSeq without dashes
            //    This is the analog of the actual reference location, but for the
            //    query.
            size_t currentRealQueryStart = currentQueryPos;
            while (tagSeq[currentRealQueryStart] == '-')
            {
                ++currentRealQueryStart;
            }
            size_t currentRealQueryPos = currentRealQueryStart;
            // Finally, variables that point to the location in the reference
            size_t currentRealTargetStart = static_cast<size_t>(referenceStart);
            size_t currentRealTargetPos = currentRealTargetStart;

            tempQueryLocation.info().clear();
            tempTargetLocation.info().clear();

            string currentType = "NONE";
            while ((currentRealQueryPos <= static_cast<size_t>(tagEnd)) &&
                    (currentRealTargetPos <= static_cast<size_t>(referenceEnd)))
            {
                // 1. check for runs of matches...
                bool foundGood = false;   // Need this so in case there are no
                // matches, we know this.
                while ((currentRealQueryPos <= static_cast<size_t>(tagEnd)) &&
                        (currentRealTargetPos <=
                         static_cast<size_t>(referenceEnd)) &&
                        (tagSeq[currentQueryPos] == refSeq[currentTargetPos]))
                {
                    foundGood = true;
                    ++currentQueryPos;
                    ++currentRealQueryPos;
                    ++currentTargetPos;
                    ++currentRealTargetPos;
                }
                // End is 1 less than currentPos;
                if (foundGood)
                {
                    anElement.type() = NoChange;
                    // Since we RC'd tagSeq and refSeq, the location is the real one.
                    tempQueryLocation.position().assign(static_cast<long>(currentRealQueryStart), static_cast<long>(currentRealQueryPos - 1));
                    tempTargetLocation.position().assign(static_cast<long>(currentRealTargetStart), static_cast<long>(currentRealTargetPos - 1));
                    tempTargetLocation.info().push_back("element_gapLength", "0");
                    tempTargetLocation.info().push_back("element_numErrors", "0");
                    anElement.queryLocation() = tempQueryLocation;
                    anElement.targetLocation() = tempTargetLocation;
                    currentFragment.push_back(anElement);
                    tempQueryLocation.info().clear();
                    tempTargetLocation.info().clear();
                    currentQueryStart = currentQueryPos;
                    // currentTargetStart = currentTargetPos;
                    currentRealQueryStart = currentRealQueryPos;
                    currentRealTargetStart = currentRealTargetPos;
                }

                // 2. check for runs of SNPs...
                foundGood = false;
                while ((currentRealQueryPos <= static_cast<size_t>(tagEnd)) &&
                        (currentRealTargetPos <= static_cast<size_t>(referenceEnd)) &&
                        (tagSeq[currentQueryPos] != '-') &&
                        (refSeq[currentTargetPos] != '-') &&
                        (tagSeq[currentQueryPos] != refSeq[currentTargetPos]))
                {
                    foundGood = true;
                    ++currentQueryPos;
                    ++currentRealQueryPos;
                    ++currentTargetPos;
                    ++currentRealTargetPos;
                }
                // End is 1 less than currentPos;
                if (foundGood)
                {
                    anElement.type() = SNPs;
                    // Since we RC'd tagSeq and refSeq, the location is the real one.
                    tempQueryLocation.position().assign(static_cast<long>(currentRealQueryStart), static_cast<long>(currentRealQueryPos - 1));
                    tempTargetLocation.position().assign(static_cast<long>(currentRealTargetStart), static_cast<long>(currentRealTargetPos - 1));
                    tempTargetLocation.info().push_back("element_gapLength", "0");
                    tempTargetLocation.info().push_back("element_numErrors", convertToString<size_t>(currentRealQueryPos - currentRealQueryStart));
                    anElement.queryLocation() = tempQueryLocation;
                    anElement.targetLocation() = tempTargetLocation;
                    currentFragment.push_back(anElement);
                    tempQueryLocation.info().clear();
                    tempTargetLocation.info().clear();
                    currentQueryStart = currentQueryPos;
                    // currentTargetStart = currentTargetPos;
                    currentRealQueryStart = currentRealQueryPos;
                    currentRealTargetStart = currentRealTargetPos;
                }

                // 3. check for runs of insertions...
                foundGood = false;
                while ((currentRealQueryPos <= static_cast<size_t>(tagEnd)) &&
                        (currentRealTargetPos <= static_cast<size_t>(referenceEnd)) &&
                        (tagSeq[currentQueryPos] != '-') &&
                        (refSeq[currentTargetPos] == '-'))
                {
                    foundGood = true;
                    ++currentQueryPos;
                    ++currentRealQueryPos;
                    ++currentTargetPos;
                }
                // End is 1 less than currentPos;
                if (foundGood)
                {
                    anElement.type() = Insertion;
                    // Since we RC'd tagSeq and refSeq, the location is the real one.
                    tempQueryLocation.position().assign(static_cast<long>(currentRealQueryStart), static_cast<long>(currentRealQueryPos - 1));
                    size_t currGap = currentRealQueryPos - currentRealQueryStart;
                    // The below is +1 on both sides since the insertions in
                    // Helicos' sms is in front of the position, whereas in
                    // MolBioLib, it happens before the position.
                    // - No +1 in front of currentRealTargetStart, since this is a
                    //   zero-length event on the target.
                    tempTargetLocation.position().assign(static_cast<long>(currentRealTargetStart), static_cast<long>(currentRealTargetPos));
                    tempTargetLocation.info().push_back("element_gapLength", convertToString<size_t>(currGap));
                    // Since we have a gap, this does not count twice towards
                    // the number of errors.
                    tempTargetLocation.info().push_back("element_numErrors", "0");
                    anElement.queryLocation() = tempQueryLocation;
                    anElement.targetLocation() = tempTargetLocation;
                    currentFragment.push_back(anElement);
                    tempQueryLocation.info().clear();
                    tempTargetLocation.info().clear();
                    currentQueryStart = currentQueryPos;
                    // currentTargetStart = currentTargetPos;
                    currentRealQueryStart = currentRealQueryPos;
                    currentRealTargetStart = currentRealTargetPos;
                }

                // 4. check for runs of deletions...
                foundGood = false;
                while ((currentRealQueryPos <= static_cast<size_t>(tagEnd)) &&
                        (currentRealTargetPos <= static_cast<size_t>(referenceEnd)) &&
                        (tagSeq[currentQueryPos] == '-') &&
                        (refSeq[currentTargetPos] != '-'))
                {
                    foundGood = true;
                    ++currentQueryPos;
                    ++currentTargetPos;
                    ++currentRealTargetPos;
                }
                // End is 1 less than currentPos;
                if (foundGood)
                {
                    anElement.type() = Deletion;
                    // Since we RC'd tagSeq and refSeq, the location is the real one.
                    // No -1 on currentRealQueryPos since this is
                    // a zero length event.
                    tempQueryLocation.position().assign(static_cast<long>(currentRealQueryStart), static_cast<long>(currentRealQueryPos));
                    size_t currGap = currentRealTargetPos - currentRealTargetStart;
                    tempTargetLocation.position().assign(static_cast<long>(currentRealTargetStart), static_cast<long>(currentRealTargetPos - 1));
                    tempTargetLocation.info().push_back("element_gapLength", convertToString<size_t>(currGap));
                    // Since we have a gap, this does not count twice towards
                    // the number of errors.
                    tempTargetLocation.info().push_back("element_numErrors", "0");
                    anElement.queryLocation() = tempQueryLocation;
                    anElement.targetLocation() = tempTargetLocation;
                    currentFragment.push_back(anElement);
                    tempQueryLocation.info().clear();
                    tempTargetLocation.info().clear();
                    currentQueryStart = currentQueryPos;
                    // currentTargetStart = currentTargetPos;
                    currentRealQueryStart = currentRealQueryPos;
                    currentRealTargetStart = currentRealTargetPos;
                }
            }
        }

        currentFragment.globalQueryLocation = tempGlobalQueryLocation;
        currentFragment.globalTargetLocation = tempGlobalTargetLocation;


        return true;

    }



    bool readPairedReads()
    {
        return false;
    }


    void setIndexFileExtensions(string queryExt = "", string targetExt = "", string extensionExt = "")
    {
        if (queryExt == "")
        {
            queryContigIndexFileExtension = "MolBioLib.HelicosSmsAlignmentFragmentReader.QueryContig.index";
        }
        else
        {
            queryContigIndexFileExtension = queryExt;
        }
        if (targetExt == "")
        {
            targetContigIndexFileExtension = "MolBioLib.HelicosSmsAlignmentFragmentReader.TargetContig.index";
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
            targetLocationIndexFileExtension = "MolBioLib.HelicosSmsAlignmentFragmentReader.TargetLocation.index";
        }
        else
        {
            targetLocationIndexFileExtension = targetLocExt;
        }
    }


    bool paired(HelicosSmsAlignmentFragmentReader& rhs)
    {
        // XXX WORK ON THIS!!! XXX
        assert(true == false);
        return false;
    }



protected:
    void createTargetLocationIndices(string filename)
    {
        // TO BE WRITTEN
        // Override this function in derived class.
        // Check that filename + extension either does not exist else is deleted.
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
        cerr << "Until HelicosSmsAlignmentFragmentReader::createTargetLocationIndices is written, HelicosSmsAlignmentFragmentReader::checkTargetLocationFile will always return 1.\n";
        return 1;
#endif
        return fileSize(tempFileName);
    }



private:
    bool readSequences, keepReads;
    Fasta *ptrSequences;

    bool isPaired;
    string tagID, referenceID, strand;
    Sequence tagSeq, refSeq;
    Sequence tagSeqWithoutDashes;
    size_t numMatch, numDel, numIns, numSub, numErrs;
    long referenceStart, referenceEnd, tagStart, tagEnd;
    double nScore;

    set<string> tagIDs;  // Set of all tagIDs.
    vector<string> theHeaderSection;

    bool hasRefHeaderSizes;
    ContigSizes refHeaderSizes;
    StringQualifiers mapChr;

    size_t headerOffset, tempIndex;
    string line;
    vector<string> tokens;

    SequenceLocation tempQueryLocation, tempTargetLocation;
    SequenceLocation tempGlobalQueryLocation, tempGlobalTargetLocation;
    SequenceLocation::position_type tempInterval;
    string tempContig;

    SequenceAlignmentFragment::alignment_element_type anElement;

    size_t currentNumErrors;

};



#endif



