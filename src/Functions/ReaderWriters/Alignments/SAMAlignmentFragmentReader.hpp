#ifndef MOLBIOLIB_SAMALIGNMENTFRAGMENTREADER_H
#define MOLBIOLIB_SAMALIGNMENTFRAGMENTREADER_H 1

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
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"
#include "src/Functions/SystemUtilities/FileModificationTime.hpp"



/** \file SAMAlignmentFragmentReader.hpp
 * The SAM format reader.
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
class SAMAlignmentFragmentReader : public SequenceAlignmentFragmentReader
{
public:

    bool aligned()
    {
        if ((RNAME == "*") || (CIGAR == "*"))
        {
            return false;
        }
        else
        {
            return true;
        }
    }


    long getMAPQ()
    {
        return MAPQ;
    }

    double alignmentScore()
    {
        return static_cast<double>(MAPQ);
    }

    long getISIZE()
    {
        return ISIZE;
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
        return true;
    }

    Fasta& getSequences()
    {
        return (*ptrSequences);
    }

    void setSequences(Fasta& theSeq)
    {
        ptrSequences = &theSeq;
    }

    Qual& getQualities()
    {
        return (*ptrQualities);
    }

    void setQualities(Qual& theQual)
    {
        ptrQualities = &theQual;
    }

    void clearSequences()
    {
        ptrSequences->clear();
    }

    void clearQualities()
    {
        ptrQualities->clear();
    }

    string getSEQ()
    {
        return SEQ;
    }
    string getCurrentSequence()
    {
        return SEQ;
    }
    string getQUAL()
    {
        return QUAL;
    }
    string getCurrentQuality()
    {
        return QUAL;
    }


    size_t numOptionalFields()
    {
        vector<string> tokens;
        splitString(optionalFields, "\t", tokens);
        return tokens.size();
    }
    string getOptionalField(size_t n)
    {
        vector<string> tokens;
        splitString(optionalFields, "\t", tokens);
        return tokens[n];
    }


    void open(string filename, string refHeaderSizesFile = "")
    {
        // We need to override this function since we also need to know if this
        // file has any header section.
        refHeaderSizes.clear();
        if (refHeaderSizesFile != "")
        {
            hasRefHeaderSizes = true;
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
        QNAMEs.clear();
        ptrQualities->clear();
        theHeaderSection.clear();
        // Find first line that does not start with @.
        while (currentFileIndex < theFileTable.size() &&
                theFileTable[currentFileIndex].substr(0, 1) == "@")
        {
            theHeaderSection.push_back(theFileTable[currentFileIndex]);
            ++currentFileIndex;
        }
        headerOffset = currentFileIndex;
        currentFileIndex = numeric_limits<size_t>::max();
        this->openIndices();
    }


    void resetToBegin(bool firstOrSecond = true,
                      bool keepReads = false)
    {
        currentFileIndex = numeric_limits<size_t>::max();
        // Assumed the header is already read in.
        if (!keepReads)
        {
            QNAMEs.clear();
            ptrSequences->clear();
            ptrQualities->clear();
        }
        doFirstOrSecond = firstOrSecond;
    }


    // firstOrSecond = true -> do first read.
    SAMAlignmentFragmentReader(bool bothFirstAndSecond = true, bool firstOrSecond = true, bool theKeepReads = true, bool theKeepQuals = true, bool theDoQueryIndex = true, bool theDoTargetIndex = true, bool theDoLocationIndex = true, bool theSequentialMode = false, float maxLoadFactor = 1.0)
    {
        hasPairedReads = false;
        keepReads = theKeepReads;
        keepQuals = theKeepQuals;
        doBothFirstAndSecond = bothFirstAndSecond;
        doFirstOrSecond = firstOrSecond;
        doQueryIndex = theDoQueryIndex;
        doTargetIndex = theDoTargetIndex;
        doLocationIndex = theDoLocationIndex;
        sequentialMode = theSequentialMode;
        hasRefHeaderSizes = false;
        pairedRead = false;
        properPair = false;
        queryUnmapped = true;
        mateUnmapped = true;
        queryStrand = false;
        mateStrand = false;
        firstRead = false;
        secondRead = false;
        notPrimary = false;
        failedQualityChecks = true;
        duplicate = true;
        string queryExt = "", targetExt = "", indexExt = "", targetLocExt = "";
        limitReads = false;
        readZN = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        this->setIndexFileExtensions(queryExt, targetExt, indexExt);
        setTargetLocationIndexFileExtension(targetLocExt);
        this->setMaxLoadFactor(maxLoadFactor);
        ptrSequences = new Fasta;
        ptrSequences->setMaxLoadFactor(maxLoadFactor);
        ptrQualities = new Qual;
        ptrQualities->setMaxLoadFactor(maxLoadFactor);
        regenerateHeaderMapsOnInsert = true;
        readSequences = true;
#ifdef DEBUG
        if (keepReads || keepQuals)
        {
            cerr << "Warning from SAMAlignmentFragmentReader!  Keeping the reads and/or the quals will make reading the SAM file much slower since it has to regenerate the header mappings on every read.  If not doing setReadSequences(false), consider batching and doing setRegenerateHeaderMapsOnInsert(false).  Continuing execution." << endl;
        }
#endif
    }

    SAMAlignmentFragmentReader(string filename, string refHeaderSizesFile = "", string queryExt = "", string targetExt = "", string indexExt = "", string targetLocExt = "", bool bothFirstAndSecond = true, bool firstOrSecond = true, bool theKeepReads = true, bool theKeepQuals = true, bool theDoQueryIndex = true, bool theDoTargetIndex = true, bool theDoLocationIndex = true, bool theSequentialMode = false, float maxLoadFactor = 1.0)
    {
        hasPairedReads = false;
        keepReads = theKeepReads;
        keepQuals = theKeepQuals;
        doBothFirstAndSecond = bothFirstAndSecond;
        doFirstOrSecond = firstOrSecond;
        doQueryIndex = theDoQueryIndex;
        doTargetIndex = theDoTargetIndex;
        doLocationIndex = theDoLocationIndex;
        sequentialMode = theSequentialMode;
        pairedRead = false;
        properPair = false;
        queryUnmapped = true;
        mateUnmapped = true;
        queryStrand = false;
        mateStrand = false;
        firstRead = false;
        secondRead = false;
        notPrimary = false;
        failedQualityChecks = true;
        duplicate = true;
        limitReads = false;
        readZN = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        this->setIndexFileExtensions(queryExt, targetExt, indexExt);
        setTargetLocationIndexFileExtension(targetLocExt);
        this->setMaxLoadFactor(maxLoadFactor);
        ptrSequences = new Fasta;
        ptrSequences->setMaxLoadFactor(maxLoadFactor);
        ptrQualities = new Qual;
        ptrQualities->setMaxLoadFactor(maxLoadFactor);
        regenerateHeaderMapsOnInsert = true;
        readSequences = true;
#ifdef DEBUG
        if (keepReads || keepQuals)
        {
            cerr << "Warning from SAMAlignmentFragmentReader!  Keeping the reads and/or the quals will make reading the SAM file much slower since it has to regenerate the header mappings on every read.  If not doing setReadSequences(false), consider batching and doing setRegenerateHeaderMapsOnInsert(false).  Continuing execution." << endl;
        }
#endif
        open(filename, refHeaderSizesFile);
    }


    // If set to true, then the ZN option will be read and stored in
    // the globalQueryLocation's info, with key ZN (of type size_t).
    void setReadZN(bool newValue)
    {
        readZN = newValue;
    }


    size_t size()
    {
        size_t numEntries = theFileTable.size() - headerOffset;
        return numEntries;
    }


    void regenerateHeaderMaps()
    {
        ptrSequences->regenerateHeaderMaps();
        ptrQualities->regenerateHeaderMaps();
    }

    /* \fn size_t numErrors()
     * The problem is that with the CIGAR format, one cannot tell how many
     * SNPs are in [0-9]+M, since M means either match or mismatch.
     */
    size_t numErrors()
    {
        return currentNumErrors;
    }

    bool readAlignmentFragment(size_t n, bool& isPaired, bool& isFirst, bool& isSecond, bool& aligns)
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

        QNAME = tokens[0];

        FLAG = convertFromString<long>(tokens[1]);
        // Below are from section 2.2.2 (<flag> field) of the SAM specifications).
        pairedRead = FLAG & 0x0001;
        isPaired = pairedRead;
        if (pairedRead)
        {
            hasPairedReads = true;
        }
        properPair = FLAG & 0x0002;
        queryUnmapped = FLAG & 0x0004;
        mateUnmapped = FLAG & 0x0008;
        queryStrand = FLAG & 0x0010;
        mateStrand = FLAG & 0x0020;
        firstRead = FLAG & 0x0040;
        isFirst = firstRead;
        secondRead = FLAG & 0x0080;
        isSecond = secondRead;
        notPrimary = FLAG & 0x0100;
        failedQualityChecks = FLAG & 0x0200;
        duplicate = FLAG & 0x0400;

        RNAME = tokens[2];
        POS = convertFromString<long>(tokens[3]);
        MAPQ = convertFromString<long>(tokens[4]);
        CIGAR = tokens[5];
        MRNM = tokens[6];
        MPOS = convertFromString<long>(tokens[7]);
        ISIZE = convertFromString<long>(tokens[8]);
        if (readSequences && keepReads)
        {
            SEQ = tokens[9];
        }
        if (readSequences && keepQuals)
        {
            QUAL = tokens[10];
        }
        size_t numTokens = tokens.size();
        optionalFields = "";
        for (size_t i = 11; i < numTokens; ++i)
        {
            optionalFields += tokens[i] + "\t";
        }
        if (readZN)
        {
            ZN = "1";
            for (size_t i = 11; i < numTokens; ++i)
            {
                if ((tokens[i][0] == 'Z') && (tokens[i][1] == 'N'))
                {
                    vector<string> znTokens;
                    splitString(tokens[i], ":", znTokens);
                    ZN = znTokens[2];
                }
            }
        }

        bool foundQname = false;
        if (readSequences && keepReads)
        {
            if (QNAMEs.find(QNAME) != QNAMEs.end())
            {
                foundQname = true;
            }
            else
            {
                if (SEQ != "*")
                {
                    if (!isPaired ||
                            doBothFirstAndSecond ||
                            (firstRead && doFirstOrSecond) ||
                            (secondRead && !doFirstOrSecond))
                    {
                        QNAMEs.insert(QNAME);
                        ptrSequences->push_back(QNAME, SEQ);
                        if (regenerateHeaderMapsOnInsert)
                        {
                            ptrSequences->regenerateHeaderMaps();
                        }
                    }
#ifdef DEBUG
                    else if (notPrimary)
                    {
                        cerr << "Warning.  Did not push_back non-primary read.\n";
                    }
#endif
                }
                else
                {
                    string blankString = "";
                    Sequence tempSeq = blankString;
                    if (!isPaired ||
                            doBothFirstAndSecond ||
                            (firstRead && doFirstOrSecond) ||
                            (secondRead && !doFirstOrSecond))
                    {
                        QNAMEs.insert(QNAME);
                        ptrSequences->push_back(QNAME, tempSeq);
                        if (regenerateHeaderMapsOnInsert)
                        {
                            ptrSequences->regenerateHeaderMaps();
                        }
                    }
#ifdef DEBUG
                    else if (notPrimary)
                    {
                        cerr << "Warning.  Did not push_back non-primary read.\n";
                    }
#endif
                }
            }
        }
        if (readSequences && keepQuals && !foundQname)
        {
            // Must check again, since keepReads may be false
            // while keepQuals is true.
            if (keepReads || (QNAMEs.find(QNAME) == QNAMEs.end()))
            {
                if (QUAL != "*")
                {
                    if (!isPaired ||
                            doBothFirstAndSecond ||
                            (firstRead && doFirstOrSecond) ||
                            (secondRead && !doFirstOrSecond))
                    {
                        QNAMEs.insert(QNAME);
                        ptrQualities->push_back(QNAME, QUAL);
                        if (regenerateHeaderMapsOnInsert)
                        {
                            ptrQualities->regenerateHeaderMaps();
                        }
                    }
#ifdef DEBUG
                    else if (notPrimary)
                    {
                        cerr << "Warning.  Did not push_back non-primary quality scores.\n";
                    }
#endif
                }
                else
                {
                    string tempQual = "";
                    if (!isPaired ||
                            doBothFirstAndSecond ||
                            (firstRead && doFirstOrSecond) ||
                            (secondRead && !doFirstOrSecond))
                    {
                        QNAMEs.insert(QNAME);
                        ptrQualities->push_back(QNAME, tempQual);
                        if (regenerateHeaderMapsOnInsert)
                        {
                            ptrQualities->regenerateHeaderMaps();
                        }
                    }
#ifdef DEBUG
                    else if (notPrimary)
                    {
                        cerr << "Warning.  Did not push_back non-primary quality scores.\n";
                    }
#endif
                }
            }
        }

        long seqLocOrigin = tempQueryLocation.origin();
        // -1 below since POS starts from 1.
        POS += seqLocOrigin - 1;
        MPOS += seqLocOrigin - 1;

        tempGlobalQueryLocation.info().clear();
        tempGlobalTargetLocation.info().clear();

        tempGlobalQueryLocation.contig() = QNAME;
        tempGlobalQueryLocation.strand() = "+";
        if (readZN)
        {
            tempGlobalQueryLocation.info().push_back("ZN", ZN);
        }
        tempGlobalTargetLocation.contig() = RNAME;
        // Yes, we want a target-centric strand, so we are flipping things.
        // SAM convention is 1 = reversed strand.
        if (queryStrand)
        {
            tempGlobalTargetLocation.strand() = "-";
        }
        else
        {
            tempGlobalTargetLocation.strand() = "+";
        }


        // Only parse if aligned.
        if (CIGAR != "*")
        {

            aligns = true;

            if (hasRefHeaderSizes)
            {
                tempGlobalTargetLocation.info().push_back("target_contig_length", convertToString<size_t>(refHeaderSizes.getSize(RNAME)));
            }
            tempGlobalTargetLocation.info().push_back("MAPQ", convertToString<long>(MAPQ));
            tempGlobalTargetLocation.info().push_back("FLAG", convertToString<long>(FLAG));
            tempGlobalTargetLocation.info().push_back("MRNM", MRNM);
            tempGlobalTargetLocation.info().push_back("ISIZE", convertToString<long>(ISIZE));
            tempGlobalTargetLocation.info().push_back("MPOS", convertToString<long>(MPOS));
            tempGlobalTargetLocation.info().push_back("optional_fields", optionalFields);

            vector< pair<size_t, char> > theOps;
            pair<size_t, char> anOp;
            size_t cigarSize = CIGAR.size();
            size_t currCigarPos = 0;
            string currSize = "";
            // Below used only in NoPoly mode.
            size_t cigarNumBases = 0;

            while (currCigarPos < cigarSize)
            {
                char& currChar = CIGAR[currCigarPos];
                if ((currChar != 'M') && (currChar != 'I') && (currChar != 'D') &&
                        (currChar != 'N') && (currChar != 'S') && (currChar != 'H') &&
                        (currChar != 'P') && (currChar != 'X') && (currChar != '='))
                {
                    currSize += currChar;
                }
                else
                {
                    // Got operation character.
                    anOp.first = convertFromString<size_t>(currSize);
                    currSize = "";   // Reset the string, so as to
                    // not keep adding on.
                    // Only count things which happen on the query contig.
                    if ((currChar == 'M') || (currChar == 'I') ||
                            (currChar == 'N') ||
                            (currChar == '=') || (currChar == 'X') )
                    {
                        cigarNumBases += anOp.first;
                    }
                    anOp.second = currChar;
                    theOps.push_back(anOp);
                }
                ++currCigarPos;
            }



            if (readerMode == NoPoly)
            {
                // Much of below copied from Normal section.
                // Please see comments therein.
                size_t queryStart = 1, queryEnd = cigarNumBases;
#ifdef PROG_DEBUG
                assert(queryEnd >= queryStart);
#endif
                tempInterval.assign(queryStart + seqLocOrigin - 1, queryEnd + seqLocOrigin - 1);
                tempGlobalQueryLocation.position() = tempInterval;
                tempGlobalQueryLocation.info().push_back("query_contig_length", convertToString<size_t>(queryEnd - queryStart + 1));
                tempInterval.assign(POS, POS + (queryEnd - queryStart));
                tempGlobalTargetLocation.position() = tempInterval;

                currentFragment.globalQueryLocation = tempGlobalQueryLocation;
                currentFragment.globalTargetLocation = tempGlobalTargetLocation;
                tempQueryLocation.info().clear();
                tempTargetLocation.info().clear();
                anElement.type() = NoChange;
                tempTargetLocation = tempGlobalTargetLocation;
                tempTargetLocation.info().push_back("element_gapLength", "0");
                anElement.queryLocation() = tempGlobalQueryLocation;
                anElement.targetLocation() = tempTargetLocation;
                currentFragment.push_back(anElement);

            }
            else
            {
                // readerMode == Normal
                // First, modify the global query start and stops based on whether
                // there are soft clips at the ends.  (Hard clips do not appear.)
                // We make queryStart = 1 instead of 0, since a -1 to queryStart
                // can occur, so this avoids the wraparound problem with size_t.
                size_t queryStart = 1, queryEnd = SEQ.size();
                if ((SEQ == "*") || !(readSequences && keepReads))
                {
                    queryEnd = cigarNumBases;
                }
                size_t numOps = theOps.size();
                bool foundNonS = false;
#ifdef DEBUG
                bool foundSH = false;
#endif
                size_t first_SH = numeric_limits<size_t>::max(),
                       last_SH = numOps;
                for (size_t i = 0; i < numOps && !foundNonS; ++i)
                {
                    if (theOps[i].second == 'S')
                    {
                        queryStart += theOps[i].first;
                    }
                    else if (theOps[i].second != 'H')
                    {
                        foundNonS = true;
                    }
                    if ((theOps[i].second == 'S') || (theOps[i].second == 'H'))
                    {
                        first_SH = i;
#ifdef DEBUG
                        foundSH = true;
#endif
                    }
                }
                if (queryStart > 1)
                {
                    tempGlobalQueryLocation.info().push_back("start_soft_clip_length", convertToString<size_t>(queryStart - 1));
                }

                if (numOps > 0)
                {
                    foundNonS = false;
                    // Must use long.  If size_t iIndex, then w/iIndex = 0, --iIndex
                    // does a wraparound for size_t of unsigned type.
                    for (long iIndex = numOps - 1; iIndex >= 0 && !foundNonS; --iIndex)
                    {
                        size_t i = static_cast<size_t>(iIndex);
                        if (theOps[i].second == 'S')
                        {
                            queryEnd -= theOps[i].first;
                        }
                        else if (theOps[i].second != 'H')
                        {
                            foundNonS = true;
                        }
                        if ((theOps[i].second == 'S') || (theOps[i].second == 'H'))
                        {
                            last_SH = i;
#ifdef DEBUG
                            foundSH = true;
#endif
                        }
                    }
                    if ((SEQ.size() - queryEnd) > 0)
                    {
                        tempGlobalQueryLocation.info().push_back("end_soft_clip_length", convertToString<size_t>(SEQ.size() - queryEnd));
                    }

#ifdef DEBUG
                    assert(queryStart <= queryEnd);

                    if (foundSH && (last_SH > (first_SH + 1)))
                    {
                        for (size_t i = first_SH + 1; i < last_SH; ++i)
                            if ((theOps[i].second == 'S') || (theOps[i].second == 'H'))
                            {
                                cerr << "Warning!  A soft or hard clip was found in the middle of an alignment, " << CIGAR << ".  These clips are ignored.  Continuing.\n";
                            }
                    }
#endif
                }

                // - 1 below since origin starts at 0.
                tempInterval.assign(queryStart + seqLocOrigin - 1, queryEnd + seqLocOrigin - 1);
                tempGlobalQueryLocation.position() = tempInterval;
                // Below used in Broad's aligner.
#ifdef PROG_DEBUG
                assert(queryEnd >= queryStart);
#endif
                tempGlobalQueryLocation.info().push_back("query_contig_length", convertToString<size_t>(queryEnd - queryStart + 1));

                // Next, modify the global target start and stops based on whether
                // there are hard or soft clips at the ends.  POS is from the examples
                // in the SAM paper the first non-clipped location in the reference.
                long targetEnd = POS - 1;
                for (size_t i = 0; i < numOps; ++i)
                {
                    char& currOp = theOps[i].second;
                    if ((currOp == 'M') || (currOp == 'D') || (currOp == 'N') ||
                            (currOp == '=') || (currOp == 'X') )
                    {
                        targetEnd += static_cast<long>(theOps[i].first);
                    }
                }
                tempInterval.assign(POS, targetEnd);
                tempGlobalTargetLocation.position() = tempInterval;

                currentFragment.globalQueryLocation = tempGlobalQueryLocation;
                currentFragment.globalTargetLocation = tempGlobalTargetLocation;

                // - 1 below since origin starts at 0.
                long currQueryPos = static_cast<long>(queryStart + seqLocOrigin - 1);
                long currTargetPos = POS;
                tempQueryLocation.contig() = tempGlobalQueryLocation.contig();
                tempQueryLocation.strand() = tempGlobalQueryLocation.strand();
                tempTargetLocation.contig() = tempGlobalTargetLocation.contig();
                tempTargetLocation.strand() = tempGlobalTargetLocation.strand();
                for (size_t i = 0; i < numOps; ++i)
                {
                    size_t& currSize = theOps[i].first;
#ifdef DEBUG
                    if (currSize == 0)
                    {
                        cerr << "Error!  Length of CIGAR operation should be > 0.  Exiting.\n";
                        assert(true == false);
                    }
#endif
                    char& currOp = theOps[i].second;
                    tempQueryLocation.info().clear();
                    tempTargetLocation.info().clear();

                    switch (currOp)
                    {
                        case 'M' :
                        case '=' :
                        case 'X' :
                            if (currOp == 'M')
                            {
                                anElement.type() = MatchOrMismatch;
                            }
                            // We do not add to currentNumErrors here because we do
                            // not know how many SNPs are in this run of M.
                            else if (currOp == '=')
                            {
                                anElement.type() = NoChange;
                            }
                            else if (currOp == 'X')
                            {
                                anElement.type() = SNPs;
                                currentNumErrors += currSize;
                            }
                            tempQueryLocation.position().assign(currQueryPos,
                                                                currQueryPos + (static_cast<long>(currSize)-1));
                            tempTargetLocation.position().assign(currTargetPos,
                                                                 currTargetPos + (static_cast<long>(currSize)-1));
                            currQueryPos += currSize;
                            currTargetPos += currSize;
                            tempTargetLocation.info().push_back("element_gapLength", "0");
                            if (currOp == 'X')
                            {
                                tempTargetLocation.info().push_back("element_numErrors", convertToString<size_t>(currSize));
                            }
                            anElement.queryLocation() = tempQueryLocation;
                            anElement.targetLocation() = tempTargetLocation;
                            currentFragment.push_back(anElement);
                            break;
                        case 'I' :
                            anElement.type() = Insertion;
                            tempQueryLocation.position().assign(currQueryPos,
                                                                currQueryPos + (static_cast<long>(currSize)-1));
                            tempTargetLocation.position().assign(currTargetPos, currTargetPos);
                            tempTargetLocation.info().push_back("element_gapLength", convertToString<size_t>(currSize));
                            tempTargetLocation.info().push_back("element_numErrors", "0");
                            currQueryPos += currSize;
                            currentNumErrors += currSize;
                            anElement.queryLocation() = tempQueryLocation;
                            anElement.targetLocation() = tempTargetLocation;
                            currentFragment.push_back(anElement);
                            break;
                        case 'D' :
                            anElement.type() = Deletion;
                            tempQueryLocation.position().assign(currQueryPos, currQueryPos);
                            tempTargetLocation.position().assign(currTargetPos,
                                                                 currTargetPos + (static_cast<long>(currSize)-1));
                            tempTargetLocation.info().push_back("element_gapLength", convertToString<size_t>(currSize));
                            tempTargetLocation.info().push_back("element_numErrors", "0");
                            currTargetPos += currSize;
                            currentNumErrors += currSize;
                            anElement.queryLocation() = tempQueryLocation;
                            anElement.targetLocation() = tempTargetLocation;
                            currentFragment.push_back(anElement);
                            break;
                        case 'N' :
                            // Skipped region on the reference
                            currTargetPos += currSize;
                            currentNumErrors += currSize;
                            break;
                        case 'P' :
                        case 'S' :
                        case 'H' :
                            if (currOp == 'P')
                            {
                                anElement.type() = Padding;
                            }
                            else if (currOp == 'S')
                            {
                                anElement.type() = SoftClip;
                            }
                            else if (currOp == 'H')
                            {
                                anElement.type() = HardClip;
                            }
                            tempQueryLocation.position().assign(currQueryPos, currQueryPos);
                            tempTargetLocation.position().assign(currTargetPos, currTargetPos);
                            tempTargetLocation.info().push_back("element_silent_length", convertToString<size_t>(currSize));
                            anElement.queryLocation() = tempQueryLocation;
                            anElement.targetLocation() = tempTargetLocation;
                            if ((i <= first_SH) || (i >= last_SH))
                            {
                                currentFragment.push_back(anElement);
                            }
                            break;
                        default :
#ifdef DEBUG
                            cerr << "This should never happen.  Bad CIGAR operation.  Exiting.\n";
                            assert(true == false);
#endif
                            break;
                    }
                }
            }

        }
        else
        {
            aligns = false;
        }

        return true;
    }


    bool readAlignmentFragment(size_t n)
    {
        bool isPaired, isFirst, isSecond, aligns;
        return readAlignmentFragment(n, isPaired, isFirst, isSecond, aligns);
    }


    bool readNextAlignmentFragment()
    {
        bool isPaired, firstResult, secondResult, aligns;
        while (true)
        {
            ++currentFileIndex;
            bool returnVal = false;
            if (!limitReads)
            {
                returnVal = readAlignmentFragment(currentFileIndex, isPaired, firstResult, secondResult, aligns);
            }
            else
            {
                bool goodRead = false;
                while (!goodRead)
                {
                    returnVal = readAlignmentFragment(currentFileIndex);
                    if (!returnVal)
                    {
                        break;
                    }
                    if (returnVal)
                    {
                        goodRead = checkGoodAlign();
                    }
                    if (!goodRead)
                    {
                        ++currentFileIndex;
                    }
                }
            }

            if (returnVal)
            {
                if ((!isPaired || doBothFirstAndSecond) && aligns)
                {
                    return true;
                }
                if (doFirstOrSecond && firstResult && aligns)
                {
                    return true;
                }
                else if (!doFirstOrSecond && secondResult && aligns)
                {
                    return true;
                }
            }
            else
            {
                return false;
            }
        }
    }


    bool readPreviousAlignmentFragment()
    {
        bool isPaired, firstResult, secondResult, aligns;
        while (true)
        {
            --currentFileIndex;
#ifdef PROG_DEBUG
            // Below - doing -- on beginning of file or first entry.
            assert((currentFileIndex != numeric_limits<size_t>::max()) &&
                   (currentFileIndex != (numeric_limits<size_t>::max()-1)));
#endif
            bool returnVal = false;
            if (!limitReads)
            {
                returnVal = readAlignmentFragment(currentFileIndex, isPaired, firstResult, secondResult, aligns);
            }
            else
            {
                bool goodRead = false;
                while (!goodRead)
                {
                    returnVal = readAlignmentFragment(currentFileIndex);
                    if (!returnVal)
                    {
                        break;
                    }
                    if (returnVal)
                    {
                        goodRead = checkGoodAlign();
                    }
                    if (!goodRead)
                    {
                        if (currentFileIndex == 0)
                        {
                            returnVal = false;
                            break;
                        }
                        else
                        {
                            --currentFileIndex;
                        }
                    }
                }
            }

            if (returnVal)
            {
                if ((!isPaired || doBothFirstAndSecond) && aligns)
                {
                    return true;
                }
                if (doFirstOrSecond && firstResult && aligns)
                {
                    return true;
                }
                else if (!doFirstOrSecond && secondResult && aligns)
                {
                    return true;
                }
            }
            else
            {
                return false;
            }
        }
    }


    bool readPairedReads()
    {
        return hasPairedReads;
    }


    void setIndexFileExtensions(string queryExt = "", string targetExt = "", string extensionExt = "")
    {
        if (queryExt == "")
        {
            queryContigIndexFileExtension = "MolBioLib.SAMAlignmentFragmentReader.QueryContig.index";
        }
        else
        {
            queryContigIndexFileExtension = queryExt;
        }
        if (targetExt == "")
        {
            targetContigIndexFileExtension = "MolBioLib.SAMAlignmentFragmentReader.TargetContig.index";
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
            targetLocationIndexFileExtension = "MolBioLib.SAMAlignmentFragmentReader.TargetLocation.index";
        }
        else
        {
            targetLocationIndexFileExtension = targetLocExt;
        }
    }


    bool paired(SAMAlignmentFragmentReader& rhs)
    {
        // XXX WORK ON THIS!!! XXX
        assert(true == false);
        return false;
    }


    bool isFirstRead()
    {
        return firstRead;
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
        cerr << "Until SAMAlignmentFragmentReader::createTargetLocationIndices is written, SAMAlignmentFragmentReader::checkTargetLocationFile will always return 1.\n";
        return 1;
#endif
        return fileSize(tempFileName);
    }



private:
    bool hasPairedReads;
    string QNAME, RNAME, CIGAR, MRNM, SEQ, QUAL, optionalFields;
    string ZN;
    set<string> QNAMEs;  // Set of all QNAMES.
    long FLAG, POS, MAPQ, MPOS, ISIZE;
    bool pairedRead, properPair, queryUnmapped, mateUnmapped, queryStrand,
         mateStrand, firstRead, secondRead, notPrimary, failedQualityChecks,
         duplicate;
    bool keepReads, keepQuals;
    bool readZN;
    size_t headerOffset, tempIndex;
    string line;
    vector<string> tokens;

    bool hasRefHeaderSizes;
    ContigSizes refHeaderSizes;

    vector<string> theHeaderSection;
    bool doFirstOrSecond;
    bool doBothFirstAndSecond;
    Fasta *ptrSequences;
    Qual *ptrQualities;

    SequenceLocation tempQueryLocation, tempTargetLocation;
    SequenceLocation tempGlobalQueryLocation, tempGlobalTargetLocation;
    SequenceLocation::position_type tempInterval;
    string tempContig;

    SequenceAlignmentFragment::alignment_element_type anElement;

    size_t currentNumErrors;

};



#endif



