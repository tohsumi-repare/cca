#ifndef MOLBIOLIB_SEQUENCEALIGNMENTFRAGMENTREADER_H
#define MOLBIOLIB_SEQUENCEALIGNMENTFRAGMENTREADER_H 1

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
#include "src/Objects/Location.hpp"
#include "src/Objects/Table.hpp"
#include "src/Objects/Fasta.hpp"
#include "src/Objects/Qual.hpp"
#include "src/Objects/ReadOnlyStringFile.hpp"
#include "src/Objects/AlignmentFragment.hpp"
#include "src/Functions/SystemUtilities/System.hpp"
#include "src/Functions/SystemUtilities/Cstdio.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"
#include "src/Functions/SystemUtilities/FileModificationTime.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"


/** \file SequenceAlignmentFragmentReader.hpp
 * Normal = usual mode, NoPoly = no polymorphisms recorded, just return as
 * perfect matches.
 * \todo Replace <code>enum AlignmentFragmentReaderMode</code> with
 * <code>enum class AlignmentFragmentReaderMode</code> when strongly typed
 * enums work.
 */
enum AlignmentFragmentReaderMode
{
    Normal, NoPoly
};


/** \file SequenceAlignmentFragmentReader.hpp
 * Defines base class #SequenceAlignmentFragmentReader.
 * Documentation on how to use #BroadQltoutAlignmentFragmentReader here.
 * Should always be using something derived from
 * #SequenceAlignmentFragmentReader, not the class itself.
 */

/** Base class to read alignment files - not meant to be used - see docs.
 * An aligment fragment is a single entry in an alignment file.  It is true
 * that sometimes, an alignment consists of just one fragment.
 * This class is not meant to be used.  Rather a class derived from this
 * should be used where readAlignmentFragment(size_t n)
 * is defined/overridden.  The expected usage of a class derived from
 * #SequenceAlignmentFragmentReader is
 * \code
 * BroadQltoutAlignmentFragmentReader myReader("test.qltout");
 *    // BroadQltoutAlignmentFragmentReader is derived from
 *      // SequenceAlignmentFragmentReader.
 *    // Alternatively, one could do
 *    // BroadQltoutAlignmentFragmentReader myReader;
 *    // myReader.open("test.qltout");
 *
 * ...
 *
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
 * if (myReader.rcQueries()) { // reverse complement the query sequence }
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
 * There is included in the #SequenceAlignmentFragmentReader a method
 * to unify all the different ways of measuring errors.
 * \code
 * double theScore = myReader.alignmentScore();
 * \endcode
 * This is highly dependent on the aligner and alignment format used.  See
 * the documentation of each type of reader to determing what this returns.
 *
 * One can ask whether a paired-end read was read by doing
 * bool readPaired = myReader.readPairedReads();
 *
 * Some derived objects, such as #BroadQltoutAlignmentFragmentReader and
 * #NCBIBlastnTableAlignmentFragmentReader allow calls to
 * \code
 * SequenceLocation myLoc;
 * ...  // myLoc is assigned to something
 * whichQueries = myReader.alignmentsRelatedToTargetLocation(myLoc);
 * \endcode
 *
 * Specific to some aligner formats, such as SAM, there is
 * \code
 * Fasta& reads = myReader.getSequences();
 * Qual& quals = myReadsergetQualities();
 * setRegenerateHeaderMapsOnInsert(bool newValue);
 * setReadSequences(bool newValue);
 * \endcode
 * The third may be called with a value of false to speed up
 * readAlignmentFragment, if getSequences() and getQualities() are not called
 * until after all the reads are done and regenerateHeaderMaps() is called
 * on each thereafter.  The latter may be called with a value of false to
 * speed up readAlignmentFragment, particularly if the sequence is not used.
 *
 * Since the format type is not always available, all derived classes must
 * provide the Boolean functions
 * \code
 * if (myReads.hasSequences()) { ... }
 * if (myReads.hasQualities()) { ... }
 * \endcode
 * so one knows if the previous methods can be called.
 *
 * On can alter the default extensions of the index filenames,
 * MolBioLib.SequenceAlignmentFragmentReader.QueryContigIndex and
 * MolBioLib.SequenceAlignmentFragmentReader.TargetContigIndex by passing in
 * the extensions (no leading period) as the last two parameters of the
 * constructor's arguments where one passes in the file name.  For the Broad
 * and Blast readers, one can also pass as the very final argument the
 * extension of the target location index file.  There is also the method
 * <code>setIndexFileExtensions</code> that accepts two strings to set the
 * file extensions.  Futhermore, there is the method
 * <code>setTargetLocationIndexFileExtension</code> that accepts a string to
 * set the target locations index file extension.  It does nothing in the
 * #SequenceAlignmentFragmentReader class, but does set it for the
 * #BroadQltoutAlignmentFragmentReader and
 * #NCBIBlastnTableAlignmentFragmentReader classes.
 *
 * There is also a parameter, maxLoadFactor that comes after all the Boolean
 * parameters, which may also be set via the function
 * <code>setMaxLoadFactor(float maxLoadFactor)</code>.)  See the documentation
 * for this in Fasta.hpp.  The default is 1.0.
 *
 * For development, the following function returns the index extension strings.
 * \code
 * string queryContigIndexExtension = "", targetContigIndexExtension = "",
 *        targetLocationIndexExtension = "", indexFileExtension = "";
 * myReader.getIndexExtensions(queryContigIndexExtension,
 *                             targetContigIndexExtension,
 *                             targetLocationIndexExtension,
 *                             indexFileExtension);
 * \endcode
 *
 * \todo Add functionality to theTargetLocationIndex.  Also, the virtual
 * function readFragment() is still a virtual function.  Someday, it
 * would be nice to make it a non-virtual function.  This is not
 * too important now, since the cost of the disk access for a single read is
 * much greater than the cost of looking up the virtual function table.
 */
class SequenceAlignmentFragmentReader
{
public:
    typedef SequenceLocation location_type;
    typedef location_type::position_type position_type;
    typedef location_type::info_type info_type;
    SequenceAlignmentFragment currentFragment;

    void setMaxLoadFactor(float maxLoadFactor)
    {
        theMaxLoadFactor = maxLoadFactor;
        theQueryContigIndex.max_load_factor(maxLoadFactor);
        theTargetContigIndex.max_load_factor(maxLoadFactor);
    }

    void setSequentialMode(bool theSequentialMode)
    {
        sequentialMode = theSequentialMode;
    }

    SequenceAlignmentFragmentReader(bool theDoQueryIndex = true,
                                    bool theDoTargetIndex = true,
                                    bool theDoLocationIndex = true,
                                    bool theSequentialMode = false,
                                    float maxLoadFactor = 1.0)
    {
        indicesOpened = false;
        senseStrandString = "+";
        antisenseStrandString = "-";
        doQueryIndex = theDoQueryIndex;
        doTargetIndex = theDoTargetIndex;
        doLocationIndex = theDoLocationIndex;
        sequentialMode = theSequentialMode;
        this->setIndexFileExtensions();
        limitReads = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        regenerateHeaderMapsOnInsert = true;
        readSequences = true;
        readerMode = Normal;
        setMaxLoadFactor(maxLoadFactor);
    }


    virtual ~SequenceAlignmentFragmentReader()
    {
        close();
        if (doQueryIndex)
        {
            theQueryContigIndex.clear();
        }
        if (doTargetIndex)
        {
            theTargetContigIndex.clear();
        }
    }


    void open(string filename)
    {
        theFileName = filename;
        theFileTable.setSequentialMode(sequentialMode);
        theFileTable.open(filename);
        currentFileIndex = numeric_limits<size_t>::max();
        openIndices();
    }

    virtual void close()
    {
        indicesOpened = false;
        limitReads = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        theFileTable.close();
        if (doQueryIndex)
        {
            ifpQueryContigIndex.close();
        }
        if (doTargetIndex)
        {
            ifpTargetContigIndex.close();
        }
    }

    SequenceAlignmentFragmentReader(string filename, string queryExt = "",
                                    string targetExt = "",
                                    string extensionExt = "",
                                    bool theDoQueryIndex = true,
                                    bool theDoTargetIndex = true,
                                    bool theDoLocationIndex = true,
                                    bool theSequentialMode = false,
                                    float maxLoadFactor = 1.0)
    {
        theFileName = filename;
        doQueryIndex = theDoQueryIndex;
        doTargetIndex = theDoTargetIndex;
        doLocationIndex = theDoLocationIndex;
        sequentialMode = theSequentialMode;
        this->setIndexFileExtensions(queryExt, targetExt, extensionExt);
        senseStrandString = "0";
        antisenseStrandString = "1";
        limitReads = false;
        maxErrors = numeric_limits<size_t>::max();
        maxAlignments = numeric_limits<size_t>::max();
        regenerateHeaderMapsOnInsert = true;
        readSequences = true;
        readerMode = Normal;
        setMaxLoadFactor(maxLoadFactor);
        open(filename);
    }


    void setReaderMode(AlignmentFragmentReaderMode newMode)
    {
        readerMode = newMode;
    }


    virtual size_t numErrors()
    {
        assert(true == false);  // Must override this function in derived class.
        return numeric_limits<size_t>::max();
    }

    // Highly aligner format dependent score.
    virtual double alignmentScore()
    {
        assert(true == false);  // Must override this function in derived class.
        return 0.0;
    }

    virtual bool readPairedReads()
    {
        return false;
    }

    virtual bool readAlignmentFragment(size_t n)
    {
        assert(true == false);  // Must override this function in derived class.
        return false;
    }


    virtual bool aligned()
    {
        // Typically want to override this in derived class, unless the alignment
        // format only consists of aligned reads, such as the Broad's qltout.
        return true;
    }

    bool paired(SequenceAlignmentFragmentReader& rhs)
    {
        assert(true == false);  // Must override this function in derived class.
        return false;
    }


    virtual bool isFirstRead()
    {
        assert(true == false);  // Must override this function in derived class.
        return false;
    }


    bool operator[](size_t n)
    {
#ifdef PROG_DEBUG
        assert(n >= 0 && n < size());
#endif
        return readAlignmentFragment(n);
    }


    virtual void filePointerResetToBegin()
    {
        theFileTable.resetToBegin();
    }


    void resetToBegin()
    {
        currentFileIndex = numeric_limits<size_t>::max();
    }

    SequenceAlignmentFragment getFragment()
    {
        return currentFragment;
    }
    SequenceAlignmentFragment& operator()()
    {
        return currentFragment;
    }

    // We really have no choice but make this virtual on this one...
    virtual size_t size()
    {
        return theFileTable.size();
    }

    size_t getIndex()
    {
        return currentFileIndex;
    }



    StringTable getQueryContigs()
    {
#ifdef DEBUG
        assert(doQueryIndex);
#endif
        StringTable result;
        for (unordered_map<string, ifstream::pos_type>::iterator i = theQueryContigIndex.begin(); i != theQueryContigIndex.end(); ++i)
        {
            result.push_back(i->first);
        }
        return result;
    }

    size_t numberAlignmentFragmentsRelatedToQueryContig(string& theContig)
    {
#ifdef DEBUG
        assert(doQueryIndex && indicesOpened);
#endif
        if (theQueryContigIndex.find(theContig) == theQueryContigIndex.end())
        {
            return 0;
        }
        ifpQueryContigIndex.clear();
        ifpQueryContigIndex.seekg(theQueryContigIndex[theContig]);
        string line;
        getline(ifpQueryContigIndex, line);
        return convertFromString<size_t>(line);
    }

    IndexTable alignmentFragmentsRelatedToQueryContig(string& theContig)
    {
#ifdef DEBUG
        assert(doQueryIndex && indicesOpened);
#endif
        IndexTable theResult;
        size_t numAligns = numberAlignmentFragmentsRelatedToQueryContig(theContig);
        string line;
        for (size_t i = 0; i < numAligns; ++i)
        {
            getline(ifpQueryContigIndex, line);
            theResult.push_back(convertFromString<size_t>(line));
        }

        return theResult;
    }


    StringTable getTargetContigs()
    {
#ifdef DEBUG
        assert(doTargetIndex);
#endif
        StringTable result;
        for (unordered_map<string, ifstream::pos_type>::iterator i = theTargetContigIndex.begin(); i != theTargetContigIndex.end(); ++i)
        {
            result.push_back(i->first);
        }
        return result;
    }

    size_t numberAlignmentFragmentsRelatedToTargetContig(string& theContig)
    {
#ifdef DEBUG
        assert(doTargetIndex && indicesOpened);
#endif
        if (theTargetContigIndex.find(theContig) == theTargetContigIndex.end())
        {
            return 0;
        }
        ifpTargetContigIndex.clear();
        ifpTargetContigIndex.seekg(theTargetContigIndex[theContig]);
        string line;
        getline(ifpTargetContigIndex, line);
        return convertFromString<size_t>(line);
    }

    IndexTable alignmentFragmentsRelatedToTargetContig(string& theContig)
    {
#ifdef DEBUG
        assert(doTargetIndex && indicesOpened);
#endif
        IndexTable theResult;
        size_t numAligns = numberAlignmentFragmentsRelatedToTargetContig(theContig);
        string line;
        for (size_t i = 0; i < numAligns; ++i)
        {
            getline(ifpQueryContigIndex, line);
            theResult.push_back(convertFromString<size_t>(line));
        }

        return theResult;
    }


    virtual size_t numberAlignmentFragmentsRelatedToTargetLocation(location_type& theLoc)
    {
        assert(true == false);   // Must override this function in derived class.
        // Only specific alignment readers know the type of
        // location to be used.  Some may not, for example, be
        // interval-based but rather position-based, in which
        // case this function has to be treated differently.
    }
    virtual IndexTable alignmentsRelatedToTargetLocation(location_type& theLoc)
    {
        assert(true == false);   // Must override this function in derived class.
        // Only specific alignment readers know the type of
        // location to be used.  Some may not, for example, be
        // interval-based but rather position-based, in which
        // case this function has to be treated differently.
    }

    virtual bool rcQueries()
    {
        assert(true == false);  // Must override this function in derived class.
    }

    virtual bool hasSequences()
    {
        assert(true == false);  // Must override this function in derived class.
    }

    virtual bool hasQualities()
    {
        assert(true == false);  // Must override this function in derived class.
    }

    virtual Fasta& getSequences()
    {
        assert(true == false);  // Must override this function in derived class.
    }

    virtual void setSequences(Fasta& theSequences)
    {
        assert(true == false);  // Must override this function in derived class.
    }

    virtual Qual& getQualities()
    {
        assert(true == false);  // Must override this function in derived class.
    }

    virtual void setQualities(Qual& theQualities)
    {
        assert(true == false);  // Must override this function in derived class.
    }

    virtual void clearSequences()
    {
        // Override in derived class, if needed.
    }

    virtual void clearQualities()
    {
        // Override in derived class, if needed.
    }

    void setRegenerateHeaderMapsOnInsert(bool newValue)
    {
        regenerateHeaderMapsOnInsert = newValue;
    }

    void setReadSequences(bool newValue)
    {
        readSequences = newValue;
    }

    virtual void setReadZN(bool newValue)
    {
        // Do nothing.  Override in derived class, if needed.
    }

    virtual void regenerateHeaderMaps()
    {
        // Override this in derived class, if needed.
        return;
    }

    virtual string getCurrentSequence()
    {
        assert(true == false);  // Must override this function in derived class.
    }

    virtual string getCurrentQuality()
    {
        assert(true == false);  // Must override this function in derived class.
    }

    void setMaxErrors(size_t maxErrs = numeric_limits<size_t>::max())
    {
        if (maxErrs == numeric_limits<size_t>::max())
        {
            limitReads = false;
        }
        else
        {
            limitReads = true;
        }
        maxErrors = maxErrs;
    }

    void setMaxAlignments(size_t maxAligns = numeric_limits<size_t>::max())
    {
        if (maxAligns == numeric_limits<size_t>::max())
        {
            limitReads = false;
        }
        else
        {
            limitReads = true;
        }
        maxAlignments = maxAligns;
    }

    bool checkGoodAlign()
    {
        if (!limitReads)
        {
            return true;
        }
        bool result = true;
        if ((maxErrors != numeric_limits<size_t>::max()) &&
                (numErrors() > maxErrors))
        {
            result = false;
        }
        if (maxAlignments != numeric_limits<size_t>::max())
        {
            string currContig = currentFragment.globalQueryLocation.contig();
            if (numberAlignmentFragmentsRelatedToQueryContig(currContig) >
                    maxAlignments)
            {
                result = false;
            }
        }
        return result;
    }

    virtual bool readNextAlignmentFragment()
    {
        ++currentFileIndex;
        if (!limitReads)
        {
            return readAlignmentFragment(currentFileIndex);
        }
        else
        {
            bool goodRead = false;
            while (!goodRead)
            {
                if (!readAlignmentFragment(currentFileIndex))
                {
                    return false;
                }
                else
                {
                    goodRead = checkGoodAlign();
                }
                if (!goodRead)
                {
                    ++currentFileIndex;
                }
            }
            return true;
        }
    }

    bool operator++()
    {
        return readNextAlignmentFragment();
    }

    bool operator++(int)
    {
        return readNextAlignmentFragment();
    }

    virtual bool readPreviousAlignmentFragment()
    {
        --currentFileIndex;
#ifdef PROG_DEBUG
        // Below - doing -- on beginning of file or first entry.
        assert((currentFileIndex != numeric_limits<size_t>::max()) &&
               (currentFileIndex != (numeric_limits<size_t>::max()-1)));
#endif
        if (!limitReads)
        {
            return readAlignmentFragment(currentFileIndex);
        }
        else
        {
            bool goodRead = false;
            while (!goodRead)
            {
                if (!readAlignmentFragment(currentFileIndex))
                {
                    return false;
                }
                else
                {
                    goodRead = checkGoodAlign();
                }
                if (!goodRead)
                {
                    // Below means we are at the beginning and so found no good reads
                    if (currentFileIndex == 0)
                    {
                        return false;
                    }
                    --currentFileIndex;
                }
            }
            return true;
        }
    }

    bool operator--()
    {
        return readPreviousAlignmentFragment();
    }

    bool operator--(int)
    {
        return readPreviousAlignmentFragment();
    }




    virtual void setIndexFileExtensions(string queryExt = "", string targetExt = "", string extensionExt = "")
    {
        if (queryExt == "")
        {
            queryContigIndexFileExtension = "MolBioLib.SequenceAlignmentFragmentReader.QueryContig.index";
        }
        else
        {
            queryContigIndexFileExtension = queryExt;
        }
        if (targetExt == "")
        {
            targetContigIndexFileExtension = "MolBioLib.SequenceAlignmentFragmentReader.TargetContig.index";
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



    virtual void setTargetLocationIndexFileExtension(string targetLocExt = "")
    {
        if (targetLocExt == "")
        {
            targetLocationIndexFileExtension = "MolBioLib.SequenceAlignmentFragmentReader.TargetLocation.index";
        }
        else
        {
            targetLocationIndexFileExtension = targetLocExt;
        }
    }


    void getIndexExtensions(string& queryContigExt, string& targetContigExt,
                            string& targetLocExt, string& indexFileExt)
    {
        queryContigExt = queryContigIndexFileExtension;
        targetContigExt = targetContigIndexFileExtension;
        targetLocExt = targetLocationIndexFileExtension;
        indexFileExt = indexFileExtension;
    }

protected:
    bool limitReads;
    size_t maxErrors, maxAlignments;
    string theFileName;
    bool doQueryIndex, doTargetIndex, doLocationIndex;
    string queryContigIndexFileExtension, targetContigIndexFileExtension,
           targetLocationIndexFileExtension, indexFileExtension;
    string senseStrandString, antisenseStrandString;
    bool sequentialMode;
    ReadOnlyStringFile theFileTable;
    size_t currentFileIndex;
    bool regenerateHeaderMapsOnInsert;
    // Below is to stop anything derived from this from rereading sequences in
    // the alignment file.  Set to false to stop.
    bool readSequences;
    float theMaxLoadFactor;

    AlignmentFragmentReaderMode readerMode;

    unordered_map<string, ifstream::pos_type> theQueryContigIndex, theTargetContigIndex;
    Ifstream ifpQueryContigIndex, ifpTargetContigIndex;



    virtual void createTargetLocationIndices(string filename)
    {
        // Override this function in derived class.
    }


    virtual void openTargetLocationIndices(string filename)
    {
        // Override this function in derived class.
    }


    virtual ifstream::pos_type checkTargetLocationFile(string filename)
    {
        // Override this function in derived class.
        return 1;
    }



    virtual void createIndices(string filename)
    {
        AlignmentFragmentReaderMode orgMode = readerMode;
        // To get the target and query contig indices, do not need to read
        // polymorphisms.  However, the readerMode is reverted before doing the
        // positional indices, since they may depend on the polymorphisms.
        // - Note that the position indices have yet to be implemented.
        readerMode = NoPoly;

        size_t numPolys = size();

        string queryContig, targetContig;
        // pos_type queryPos, targetPos;

        unordered_map< string, set<size_t> > mapQueryContigIndex, mapTargetContigIndex;
        mapQueryContigIndex.max_load_factor(theMaxLoadFactor);
        mapTargetContigIndex.max_load_factor(theMaxLoadFactor);

        // We do not want to regenerate the header maps of the sequences and
        // qualities right now.
        bool currRegen = regenerateHeaderMapsOnInsert;
        regenerateHeaderMapsOnInsert = false;
        bool currReadSeq = readSequences;
        readSequences = false;
        if (sequentialMode)
        {
            this->filePointerResetToBegin();
        }
        for (size_t i = 0; i < numPolys; ++i)
        {
            readAlignmentFragment(i);
            if (!(this->aligned()))
            {
                continue;
            }
            if (doQueryIndex)
            {
                queryContig = currentFragment.getQueryLocation(0).contig();
                mapQueryContigIndex[queryContig].insert(i);
            }
            if (doTargetIndex)
            {
                targetContig = currentFragment.getTargetLocation(0).contig();
                mapTargetContigIndex[targetContig].insert(i);
            }

            // For interval later...
            // if (doLocationIndex) {
            //    queryPos = currentFragment.getQueryLocation().getPos();
            //    targetPos = currentFragment.getQueryLocation().getPos();
            //    ...
            // }
        }
        if (sequentialMode)
        {
            this->filePointerResetToBegin();
        }
        regenerateHeaderMapsOnInsert = currRegen;
        if (regenerateHeaderMapsOnInsert)
        {
            this->regenerateHeaderMaps();
        }
        readSequences = currReadSeq;

        Ofstream fout, foutIndex;

        if (doQueryIndex)
        {
            string tempFilename = filename + "." + queryContigIndexFileExtension;
            fout.open(tempFilename);
            string tempIndexFilename = filename + "." + queryContigIndexFileExtension + "." + indexFileExtension;
            foutIndex.open(tempIndexFilename);
            for (unordered_map< string, set<size_t> >::iterator currIndex = mapQueryContigIndex.begin(); currIndex != mapQueryContigIndex.end(); ++currIndex)
            {
                foutIndex << currIndex->first << "\t" << fout.tellp() << endl;
                set<size_t>& currentAlignSet = currIndex->second;
                fout << currentAlignSet.size() << endl;
                for (set<size_t>::iterator k = currentAlignSet.begin(); k != currentAlignSet.end(); ++k)
                {
                    fout << *k << endl;
                }
            }
            fout.close();
            foutIndex.close();
            mapQueryContigIndex.clear();
        }

        if (doTargetIndex)
        {
            string tempFilename = filename + "." + targetContigIndexFileExtension;
            fout.open(tempFilename);
            string tempIndexFilename = filename + "." + targetContigIndexFileExtension + "." + indexFileExtension;
            foutIndex.open(tempIndexFilename);
            for (unordered_map< string, set<size_t> >::iterator currIndex = mapTargetContigIndex.begin(); currIndex != mapTargetContigIndex.end(); ++currIndex)
            {
                foutIndex << currIndex->first << "\t" << fout.tellp() << endl;
                set<size_t>& currentAlignSet = currIndex->second;
                fout << currentAlignSet.size() << endl;
                for (set<size_t>::iterator k = currentAlignSet.begin(); k != currentAlignSet.end(); ++k)
                {
                    fout << *k << endl;
                }
            }
            fout.close();
            foutIndex.close();
            mapTargetContigIndex.clear();
        }
        readIndices(filename);
        readerMode = orgMode;

        if (doLocationIndex)
        {
            createTargetLocationIndices(filename);
        }
    }



    void readIndices(string filename)
    {
        if (doQueryIndex)
        {
            theQueryContigIndex.clear();
        }
        if (doTargetIndex)
        {
            theTargetContigIndex.clear();
        }

        Ifstream fin;

        string tempContig;

        string line;
        vector<string> tokens;

        if (doQueryIndex)
        {
            string tempFileName = filename + "." + queryContigIndexFileExtension + "." + indexFileExtension;
#ifdef DEBUG
            ifstream::pos_type fileCheck = fileSize(tempFileName);
            if (fileCheck == static_cast<ifstream::pos_type>(-1))
            {
                cerr << "Error in SequenceAlignmentFragmentReader::readIndices!  " << tempFileName << " could not be opened.  Exiting." << endl;
                assert(true == false);
            }
#endif
            fin.open(tempFileName);
            while (!fin.fail())
            {
                getline(fin, line);
                if (fin.fail())
                {
                    break;
                }
                splitString(line, "\t", tokens);
                tempContig = tokens[0];
                theQueryContigIndex[tempContig] = convertFromString<ifstream::pos_type>(tokens[1]);
            }
            fin.close();
            tempFileName = filename + "." + queryContigIndexFileExtension;
#ifdef DEBUG
            fileCheck = fileSize(tempFileName);
            if (fileCheck == static_cast<ifstream::pos_type>(-1))
            {
                cerr << "Error in SequenceAlignmentFragmentReader::readIndices!  " << tempFileName << " could not be opened.  Exiting." << endl;
                assert(true == false);
            }
#endif
            ifpQueryContigIndex.open(tempFileName);
        }

        if (doTargetIndex)
        {
            string tempFileName = filename + "." + targetContigIndexFileExtension + "." + indexFileExtension;
#ifdef DEBUG
            ifstream::pos_type fileCheck = fileSize(tempFileName);
            if (fileCheck == static_cast<ifstream::pos_type>(-1))
            {
                cerr << "Error in SequenceAlignmentFragmentReader::readIndices!  " << tempFileName << " could not be opened.  Exiting." << endl;
                assert(true == false);
            }
#endif
            fin.open(tempFileName);
            while (!fin.fail())
            {
                getline(fin, line);
                if (fin.fail())
                {
                    break;
                }
                splitString(line, "\t", tokens);
                tempContig = tokens[0];
                theTargetContigIndex[tempContig] = convertFromString<ifstream::pos_type>(tokens[1]);
            }
            fin.close();
            tempFileName = filename + "." + targetContigIndexFileExtension;
#ifdef DEBUG
            fileCheck = fileSize(tempFileName);
            if (fileCheck == static_cast<ifstream::pos_type>(-1))
            {
                cerr << "Error in SequenceAlignmentFragmentReader::readIndices!  " << tempFileName << " could not be opened.  Exiting." << endl;
                assert(true == false);
            }
#endif
            ifpTargetContigIndex.open(tempFileName);
        }
    }

    bool indicesOpened;

    void openIndices()
    {
        // Below is not a good test since in some cases, theFileTable is not used.
        // assert(theFileTable.isOpen() == true);
        if (indicesOpened)
        {
            cerr << "Warning openIndices() already called before.\n";
            return;
        }
        string tempFileName = theFileName + "." + queryContigIndexFileExtension;
        bool toCreate = false;
        ifstream::pos_type fileCheck = fileSize(tempFileName);
        if (doQueryIndex && (fileCheck == static_cast<ifstream::pos_type>(-1)))
        {
            toCreate = true;
        }
        else if (doQueryIndex && (fileCheck != static_cast<ifstream::pos_type>(-1))
                 && isFileOlderThan(tempFileName, theFileName))
        {
            remove(tempFileName);
            toCreate = true;
        }
        tempFileName = theFileName + "." + targetContigIndexFileExtension;
        fileCheck = fileSize(tempFileName);
        if (doTargetIndex && (fileCheck == static_cast<ifstream::pos_type>(-1)))
        {
            toCreate = true;
        }
        else if (doTargetIndex && (fileCheck != static_cast<ifstream::pos_type>(-1))
                 && isFileOlderThan(tempFileName, theFileName))
        {
            remove(tempFileName);
            toCreate = true;
        }
        tempFileName = theFileName + "." + targetLocationIndexFileExtension;
#ifdef DEBUG
        cerr << "Warning!  In SequenceAlignmentFragmentReader::openIndices(), skipping check for " << tempFileName << " since the targetLocation index routines have not been written for any of the alignment formats." << endl;
#endif
        // fileCheck = checkTargetLocationFile(tempFileName);
        // if (doLocationIndex && (fileCheck == static_cast<ifstream::pos_type>(-1)))
        //     toCreate = true;
        // XXXXXX UNCOMMENT BELOW WHEN TARGET LOCATION INDEX
        //        ROUTINES HAVE BEEN WRITTEN XXXXXX
        // else if (doLocationIndex && (fileCheck != static_cast<ifstream::pos_type>(-1))
        //                           && isFileOlderThan(tempFileName, theFileName))
        //      toCreate = true;
        if (toCreate)
        {
            createIndices(theFileName);
        }
        else
        {
            readIndices(theFileName);
        }

        // Regardless of creation, open TargetLocationIndicesFile.
        openTargetLocationIndices(theFileName);


        indicesOpened = true;
    }



};

#endif

