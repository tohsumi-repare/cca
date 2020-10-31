#ifndef MOLBIOLIB_FILESEQUENCEALIGNMENTFRAGMENTREADER_H
#define MOLBIOLIB_FILESEQUENCEALIGNMENTFRAGMENTREADER_H 1

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


#include "src/Objects/Fasta.hpp"
#include "src/Functions/ReaderWriters/Alignments/AlignmentFileTypes.hpp"
#include "src/Functions/ReaderWriters/Alignments/SequenceAlignmentFragmentReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/BroadQltoutAlignmentFragmentReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/HelicosSmsAlignmentFragmentReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/SAMAlignmentFragmentReader.hpp"



/** \file FileSequenceAlignmentFragmentReader.hpp
 * Defines class #FileSequenceAlignmentFragmentReader.
 * A single interface for the alignment readers.  This is based on the
 * Bridge design pattern of the book by Gamma et. al., "Design Patterns".
 * The implementation are the specific aligner's reader.
 */

/** Class to read alignment files.
 * An aligment fragment is a single entry in an alignment file.
 * \code
 * Fasta theReads("reads.fasta");
 * FileSequenceAlignmentFragmentReader myReader("test.qltout", "ref.headerSizes", theReads);
 * ++myReader;   // Read next polymorphism.
 *               // NOTE: THE ++'s can go on either side.
 *   // Alternatively, do myReader.readNextPolymorphism();
 * --myReader;   // Read previous polymorphism.
 *               // NOTE: THE --'s can go on either side.
 *   // Alternatively, do myReader.readPreviousPolymorphism();
 * myReader.resetToBegin();
 * // One can go back to the start of the alignments by doing the above.
 * myReader[5];   // Read sixth polymorphism.
 *    // To be more precise, one may wish to check if an error has been thrown
 *    // when DEBUG has been defined, e.g. end of file, etc. in which case,
 *    // just check the return value, i.e.
 * if (!myReader[5]) { // Treat end of file condition here. }
 *    // Can also do above with instead of myReader[5],
 *    // ++myReader and --myReader
 * // If one uses ++ or -- or

 * // Below gets the data via reference (so can be on the left side of =).
 * SequenceAlignmentFragment currFrag = myReader();
 *    // An alternative to myReader() is myReader.currentFragment
 *    // which gives direct access to the currentFragment or
 *    // use the function call myReader.getFragment() [return by value so
 *    // cannot be on left side of =).
 * \endcode
 * One can know the number of alignments and which entry is being read by
 * \code
 * size_t theSize = myReader.size();
 * size_t theEntry = myReader.getIndex();
 * \endcode
 * <b>Note that the indexing of queries assume that a polymorphism event
 * is limited to the same contig.</b>
 *
 * One can look up which queries are associated with a particular query or
 * target contig by first doing (assuming an alignment file is open):
 * \code
 * IndexTable whichQueries;
 * string n;
 * myReader.openIndices();
 * \endcode
 * followed later in the code by
 * \code
 * size_t numQueries = myReader.numberAlignmentFragmentsRelatedToQueryContig(n);
 * whichQueries = myReader.alignmentFragmentsRelatedToQueryContig(n);
 * numQueries = myReader.numberAlignmentFragmentsRelatedToTargetContig(n);
 * whichQueries = myReader.alignmentFragmentsRelatedToTargetContig(n);
 * \endcode
 * where n is set to some contig.  For the #BroadQltoutAlignmentFragmentReader,
 * n is a unsigned long where for the #NCBIBlastnTableAlignmentFragmentReader,
 * n is a string.
 *
 * Sometimes, one will wish to know whether queries found in the fasta file (or
 * in associated sequences are to be reverse-complemented against the
 * reference, for example in coverage, when the alignment is on the minus
 * strand.  To determine this, do
 * \code
 * if (myReader.rcQueries()) { // reverse complement the read }
 * \endcode
 *
 * One may wish to know all of the contigs represented in the alignments.  To
 * do so, do
 * \code
 * StringTable theQueryContigs = myReader.getQueryContigs();
 * // or
 * StringTable theTargetContigs = myReader.getTargetContigs();
 * \endcode
 * depending on which is needed.  Each contig is represented only once and
 * is in alphabetical order.
 *
 * One can ask whether a paired-end read was read by doing
 * bool readPaired = myReader.readPairedReads();
 *
 * On can alter the default extensions of the index filenames,
 * MolBioLib.SequenceAlignmentFragmentReader.QueryContigIndex and
 * MolBioLib.SequenceAlignmentFragmentReader.TargetContigIndex by passing in
 * the extensions (no leading period) as the last two parameters of the
 * constructor's arguments where one passes in the file name.
 * One can also pass as the very final argument the
 * extension of the target location index file.  There is also the method
 * <code>setIndexFileExtensions</code> that accepts two strings to set the
 * file extensions.  Futhermore, there is the method
 * <code>setTargetLocationIndexFileExtension</code> that accepts a string to
 * set the target locations index file extension.
 *
 * There is also a parameter, maxLoadFactor that comes after all of the other
 * parameters, which may also be set via the function
 * <code>setMaxLoadFactor(float maxLoadFactor)</code>.)  See the documentation
 * for this in Fasta.hpp.  The default is 1.0.
 *
 * \todo Add functionality to theTargetLocationIndex.
 */
class FileSequenceAlignmentFragmentReader
{
public:
    typedef SequenceLocation location_type;
    typedef location_type::position_type position_type;
    typedef location_type::info_type info_type;
    SequenceAlignmentFragment currentFragment;

    FileSequenceAlignmentFragmentReader(string filename,
                                        string headerSizesFilename,
                                        Fasta& reads,
                                        AlignmentFileType theAlignmentType = Null,
                                        string queryExt = "",
                                        string targetExt = "",
                                        string indexExt = "",
                                        string targetLocExt = "",
                                        bool doQueryIndex = true,
                                        bool doTargetIndex = true,
                                        bool doLocationIndex = true,
                                        bool samDoBothFirstAndSecond = true,
                                        bool samDoFirstOrSecond = true,
                                        bool keepReads = true,
                                        bool keepQuals = true,
                                        size_t maxLengthToCheck = 10000,
                                        bool sequentialMode = false,
                                        bool maxLoadFactor = 1.0,
                                        AlignmentFragmentReaderMode theMode = Normal)
    {

        if (theAlignmentType == Null)
            alignmentType = determineAlignmentFileType(filename, false,
                            maxLengthToCheck);
        else
        {
            alignmentType = theAlignmentType;
        }

        switch (alignmentType)
        {
            case Qltout:
                createBroadMapStringQualifiersFromHeaderSizes(headerSizesFilename, mapBroadChr);
                createBroadMapStringQualifiersFromFasta(reads, mapBroadQueries);
                ptrReader = new BroadQltoutAlignmentFragmentReader(filename, BroadContigStringQualifiers, mapBroadQueries, BroadContigStringQualifiers, mapBroadChr, queryExt, targetExt, indexExt, targetLocExt, doQueryIndex, doTargetIndex, doLocationIndex, sequentialMode, maxLoadFactor);
                break;
            case SAM:
                ptrReader = new SAMAlignmentFragmentReader(filename, headerSizesFilename, queryExt, targetExt, indexExt, targetLocExt, samDoBothFirstAndSecond, samDoFirstOrSecond, keepReads, keepQuals, doQueryIndex, doTargetIndex, doLocationIndex, sequentialMode, maxLoadFactor);
                ptrReader->setSequences(reads);
                ptrReader->regenerateHeaderMaps();
                break;
            case Sms:
                ptrReader = new HelicosSmsAlignmentFragmentReader(filename, headerSizesFilename, queryExt, targetExt, indexExt, targetLocExt, keepReads, doQueryIndex, doTargetIndex, doLocationIndex, sequentialMode, maxLoadFactor);
                ptrReader->setSequences(reads);
                ptrReader->regenerateHeaderMaps();
                break;
            case Null:
            default:
                cerr << "Error in FileSequenceAlignmentFragmentReader.  The alignmentType should not be Null or of unknown type.  Exiting." << endl;
                assert(true == false);
                break;
        }
        ptrReader->setReaderMode(theMode);
    }

    ~FileSequenceAlignmentFragmentReader()
    {
        ptrReader->close();
        BroadQltoutAlignmentFragmentReader *ptrBQAFR;
        HelicosSmsAlignmentFragmentReader *ptrHSAFR;
        SAMAlignmentFragmentReader *ptrSAFR;
        switch (alignmentType)
        {
            case Qltout:
                ptrBQAFR = static_cast<BroadQltoutAlignmentFragmentReader*>(ptrReader);
                delete ptrBQAFR;
                break;
            case SAM:
                ptrSAFR = static_cast<SAMAlignmentFragmentReader*>(ptrReader);
                delete ptrSAFR;
                break;
            case Sms:
                ptrHSAFR = static_cast<HelicosSmsAlignmentFragmentReader*>(ptrReader);
                delete ptrHSAFR;
                break;
            case Null:
            default:
                break;
        }
    }


    void setMaxLoadFactor(float maxLoadFactor)
    {
        ptrReader->setMaxLoadFactor(maxLoadFactor);
    }


    void setReaderMode(AlignmentFragmentReaderMode newMode)
    {
        ptrReader->setReaderMode(newMode);
    }


    void open(string filename)
    {
        ptrReader->open(filename);
    }

    void close()
    {
        ptrReader->close();
    }

    void regenerateHeaderMaps()
    {
        ptrReader->regenerateHeaderMaps();
    }

    size_t numErrors()
    {
        return ptrReader->numErrors();
    }

    bool aligned()
    {
        return ptrReader->aligned();
    }

    double alignmentScore()
    {
        return ptrReader->alignmentScore();
    }

    bool readPairedReads()
    {
        return ptrReader->readPairedReads();
    }

    bool readAlignmentFragment(size_t n)
    {
        return ptrReader->readAlignmentFragment(n);
    }


    bool paired(SequenceAlignmentFragmentReader& rhs)
    {
        return ptrReader->paired(rhs);
    }

    bool isFirstRead()
    {
        return ptrReader->isFirstRead();
    }

    bool operator[](size_t n)
    {
        return ptrReader->readAlignmentFragment(n);
    }

    void filePointerResetToBegin()
    {
        ptrReader->filePointerResetToBegin();
    }

    void resetToBegin()
    {
        ptrReader->resetToBegin();
    }

    SequenceAlignmentFragment getFragment()
    {
        return ptrReader->currentFragment;
    }
    SequenceAlignmentFragment& operator()()
    {
        return ptrReader->currentFragment;
    }

    size_t size()
    {
        return ptrReader->size();
    }

    size_t getIndex()
    {
        return ptrReader->getIndex();
    }

    bool rcQueries()
    {
        return ptrReader->rcQueries();

    }

    StringTable getQueryContigs()
    {
        return ptrReader->getQueryContigs();
    }

    size_t numberAlignmentFragmentsRelatedToQueryContig(string& theContig)
    {
        return ptrReader->numberAlignmentFragmentsRelatedToQueryContig(theContig);
    }

    IndexTable alignmentFragmentsRelatedToQueryContig(string& theContig)
    {
        return ptrReader->alignmentFragmentsRelatedToQueryContig(theContig);
    }


    StringTable getTargetContigs()
    {
        return ptrReader->getTargetContigs();
    }

    size_t numberAlignmentFragmentsRelatedToTargetContig(string& theContig)
    {
        return ptrReader->numberAlignmentFragmentsRelatedToTargetContig(theContig);
    }

    IndexTable alignmentFragmentsRelatedToTargetContig(string& theContig)
    {
        return ptrReader->alignmentFragmentsRelatedToTargetContig(theContig);
    }


    size_t numberAlignmentFragmentsRelatedToTargetLocation(location_type& theLoc)
    {
        return ptrReader->numberAlignmentFragmentsRelatedToTargetLocation(theLoc);
    }
    IndexTable alignmentsRelatedToTargetLocation(location_type& theLoc)
    {
        return ptrReader->alignmentsRelatedToTargetLocation(theLoc);
    }

    void setReadZN(bool newValue)
    {
        ptrReader->setReadZN(newValue);
    }

    void setMaxErrors(size_t maxErrs = numeric_limits<size_t>::max())
    {
        return ptrReader->setMaxErrors(maxErrs);
    }

    void setMaxAlignments(size_t maxAligns = numeric_limits<size_t>::max())
    {
        return ptrReader->setMaxAlignments(maxAligns);
    }

    bool checkGoodAlign()
    {
        return ptrReader->checkGoodAlign();
    }

    bool readNextAlignmentFragment()
    {
        return ptrReader->readNextAlignmentFragment();
    }

    bool operator++()
    {
        return ptrReader->readNextAlignmentFragment();
    }

    bool operator++(int)
    {
        return ptrReader->readNextAlignmentFragment();
    }

    bool readPreviousAlignmentFragment()
    {
        return ptrReader->readPreviousAlignmentFragment();
    }

    bool operator--()
    {
        return ptrReader->readPreviousAlignmentFragment();
    }

    bool operator--(int)
    {
        return ptrReader->readPreviousAlignmentFragment();
    }

    bool hasSequences()
    {
        return ptrReader->hasSequences();
    }

    bool hasQualities()
    {
        return ptrReader->hasQualities();
    }

    Fasta& getSequences()
    {
        return ptrReader->getSequences();
    }
    void setSequences(Fasta& theSequences)
    {
        ptrReader->setSequences(theSequences);
    }
    string getCurrentSequence()
    {
        return ptrReader->getCurrentSequence();
    }

    Qual& getQualities()
    {
        return ptrReader->getQualities();
    }
    void setQualities(Qual& theQualities)
    {
        ptrReader->setQualities(theQualities);
    }
    string getCurrentQuality()
    {
        return ptrReader->getCurrentQuality();
    }

    void clearSequences()
    {
        ptrReader->clearSequences();
    }

    void clearQualities()
    {
        ptrReader->clearQualities();
    }


    void setRegenerateHeaderMapsOnInsert(bool newValue)
    {
        ptrReader->setRegenerateHeaderMapsOnInsert(newValue);
    }

    void setReadSequences(bool newValue)
    {
        ptrReader->setReadSequences(newValue);
    }

    void setIndexFileExtensions(string queryExt = "", string targetExt = "", string extensionExt = "")
    {
        return ptrReader->setIndexFileExtensions(queryExt, targetExt, extensionExt);
    }



    void setTargetLocationIndexFileExtension(string targetLocExt = "")
    {
        return ptrReader->setTargetLocationIndexFileExtension(targetLocExt);
    }


private:
    AlignmentFileType alignmentType;
    // Below needed specifically for the Broad qltout alignment reader.
    StringQualifiers mapBroadQueries, mapBroadChr;
    SequenceAlignmentFragmentReader *ptrReader;


};

#endif

