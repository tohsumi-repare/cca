#ifndef MOLBIOLIB_FILESEQUENCEALIGNMENTFRAGMENTWRITER_H
#define MOLBIOLIB_FILESEQUENCEALIGNMENTFRAGMENTWRITER_H 1

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
#include "src/Objects/AlignmentFragment.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"


/** \file FileSequenceAlignmentFragmentWriter.hpp
 * Defines class #FileSequenceAlignmentFragmentWriter.
 */

/** Class to write alignment files
 * An aligment fragment is a single entry in an alignment file.  It is true
 * that sometimes, an alignment consists of just one fragment.
 * This class is not meant to be used.  Rather a class derived from this
 * should be used where writeAlignmentFragment
 * is defined/overridden.  The expected usage of a class derived from
 * #SequenceAlignmentFragmentWriter is
 * \code
 * SequenceAlignmentFragment theFragment;  // fill in later...
 * Fasta theReads;  // fill in later..
 * Qual theQuals;  // fill in later..
 * // ...
 * // Assuming there is a ref.headerSizes, which the output from
 * // #readsHeaderSizes with the HEADER_SIZES and KEEP_ONLY_FIRST_WORD=true 
 * // options applied on the reference,
 * FileSequenceAlignmentFragmentWriter myWriter("test.qltout",
 *                                              "ref.headerSizes",
 *                                               theReads, theQuals);
 * myWriter.write(theFragment);
 * \endcode
 */
class FileSequenceAlignmentFragmentWriter
{
public:
    typedef SequenceLocation location_type;
    typedef location_type::position_type position_type;
    typedef location_type::element_type element_type;
    typedef location_type::info_type info_type;

    FileSequenceAlignmentFragmentWriter(string filename,
                                        string headerSizesFilename,
                                        Fasta& theReads,
                                        Qual& theQuals,
                                        AlignmentFileType theAlignmentType = Null,
                                        bool samKeepReads = true,
                                        bool samKeepQual = true,
                                        size_t maxLengthToCheck = 10000) : reads(theReads), quals(theQuals)
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
                createBroadMapStringQualifiersFromHeaderSizes(headerSizesFilename,
                        mapRef);
                createBroadMapStringQualifiersFromFasta(reads, mapReads);
                ofp.open(filename);
                // true below clears mapRef and mapReads.  Not needed.
                ptrWriter = new BroadQltoutAlignmentFragmentWriter(ofp, BroadContigStringQualifiers, mapReads, BroadContigStringQualifiers, mapRef, true);
                break;
            case SAM:
                ofp.open(filename);
                ptrSAMWriter = new SAMAlignmentFragmentWriter(ofp);
                ptrSAMWriter->addFasta(theReads);
                assert(quals.size() == theReads.size());
                ptrSAMWriter->addQualities(quals);
                ptrWriter = dynamic_cast<SequenceAlignmentFragmentWriter*>(ptrSAMWriter);
                break;
            case Null:
            default:
                cerr << "Error in FileSequenceAlignmentFragmentReader.  The alignmentType should not be Null or of unknown type.  Exiting." << endl;
                assert(true == false);
                break;
        }
    }


    void close()
    {
        ofp.close();
        ptrWriter->close();
    }


    ~FileSequenceAlignmentFragmentWriter()
    {
        close();
        delete ptrWriter;
    }


    void write(SequenceAlignmentFragment& theFragment)
    {
        ptrWriter->write(theFragment);
    }


private:
    Ofstream ofp;
    AlignmentFileType alignmentType;
    SAMAlignmentFragmentWriter *ptrSAMWriter;
    SequenceAlignmentFragmentWriter *ptrWriter;
    StringQualifiers mapReads, mapRef;
    Fasta& reads;
    Qual& quals;

};

#endif

