#ifndef MOLBIOLIB_SAMALIGNMENTFRAGMENTWRITER_H
#define MOLBIOLIB_SAMALIGNMENTFRAGMENTWRITER_H 1

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
#include "src/Objects/AlignmentFragment.hpp"
#include "src/Functions/ReaderWriters/Alignments/SequenceAlignmentFragmentWriter.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"



/** \file SAMAlignmentFragmentWriter.hpp
 * Writes a SAM line in a presumably open SAM file given a fragment.
 * The usage is
 * \code
 * SequenceAlignmentFragment theFragment;   // Initialize later...
 * ostream out;   // Initialize later.  Can also be of type Ofstream, etc.
 * Fasta myReads;
 * Qual myQuals;
 * // ...
 * SAMAlignmentFragmentWriter theWriter(out, myReads, myQuals, true, true, false);
 *   // where first <code>true</code> is to print the Fasta sequences in the
 *   // SAM file, the second <code>true</code> is to print the Qual string in
 *   // the SAM file, and the last <code>false</code> means to use the regular
 *   // CIGAR string (no X's or ='s).  The last three parameters are optional
 *   // with the defaults shown.
 * \endcode
 * where the last four parameters a defaults.
 *
 * <b> Warning!  The qualities must be in the same order (with the same number
 * of entries) as the reads! </b>
 */
class SAMAlignmentFragmentWriter : public SequenceAlignmentFragmentWriter
{
public:
    SAMAlignmentFragmentWriter(ostream& ofp,
                               Fasta *ptrInputReads, Qual *ptrInputQualities,
                               bool thePrintFasta = true,
                               bool thePrintQuals = true,
                               bool theUseNewCigar = false) : out(ofp)
    {
        ptrReads = ptrInputReads;
        ptrQualities = ptrInputQualities;
        printFasta = thePrintFasta;
        printQuals = thePrintQuals;
        useNewCigar = theUseNewCigar;
    }
    SAMAlignmentFragmentWriter(ostream& ofp,
                               Fasta& inputReads, Qual& inputQualities,
                               bool theUseNewCigar = false) : out(ofp)
    {
        ptrReads = &inputReads;
        ptrQualities = &inputQualities;
        printFasta = true;
        printQuals = true;
        useNewCigar = theUseNewCigar;
    }
    SAMAlignmentFragmentWriter(ostream& ofp,
                               Fasta& inputReads,
                               int defaultQuality = -256,
                               int offset = 0,
                               bool theUseNewCigar = false) : out(ofp)
    {
        ptrReads = &inputReads;
        ptrQualities = nullptr;
        printFasta = true;
        if (defaultQuality == -256)
        {
            printQuals = false;
        }
        else
        {
#ifdef DEBUG
            assert((defaultQuality >= -5) && (defaultQuality <= 93));
#endif
            ptrQualities = new Qual(inputReads, defaultQuality, offset);
            printQuals = true;
        }
        useNewCigar = theUseNewCigar;
    }
    // Makes no sense to have the below.  Omitted.
    // SAMAlignmentFragmentWriter(ostream& out,
    //                            Qual& inputQualities,
    //                            bool theUseNewCigar = false) {
    // }
    SAMAlignmentFragmentWriter(ostream& ofp,
                               bool theUseNewCigar = false) : out(ofp)
    {
        ptrReads = nullptr;
        ptrQualities = nullptr;
        printFasta = false;
        printQuals = false;
        useNewCigar = theUseNewCigar;
    }


    void addFasta(Fasta& inputReads)
    {
        ptrReads = &inputReads;
        printFasta = true;
    }


    void addQualities(Qual& inputQualities)
    {
        ptrQualities = &inputQualities;
        printQuals = true;
    }


    void write(SequenceAlignmentFragment& theFragment)
    {
        SequenceLocation& globalQueryLocation = theFragment.globalQueryLocation;
        SequenceLocation& globalTargetLocation = theFragment.globalTargetLocation;

        SequenceLocation::element_type theOrigin = globalQueryLocation.origin();

        vector<string> words;
        // Need to only print out first word of the contig - SAM specifications!
        splitString(globalQueryLocation.contig(), " ", words);
        out << words[0] << "\t";

        StringQualifiers& targetInfo = globalTargetLocation.info();
        if (targetInfo.hasKey("FLAG"))
        {
            out << targetInfo.getQualifier("FLAG") << "\t";
        }
        else
        {
            long FLAG = 0;
            // SAM convention is - = true on flag 0x0010.
            if (globalTargetLocation.strand() == "-")
            {
                FLAG |= 0x0010;    // Sets that bit to true.
            }
            out << FLAG << "\t";
        }

        // Need to only print out first word of the contig - SAM specifications!
        splitString(globalTargetLocation.contig(), " ", words);
        out << words[0] << "\t";

        element_type POS = globalTargetLocation.position().lower();
        // Below +1 is because SAM uses 1-based index.
        POS += 1 - theOrigin;
        out << POS << "\t";

        if (targetInfo.hasKey("MAPQ"))
        {
            out << targetInfo.getQualifier("MAPQ") << "\t";
        }
        else
        {
            out << "255\t";    // According to SAM specs, if MAPQ is not available,
        }
        // use 255.

        string cigarString = "";
        size_t numElements = theFragment.size();
        element_type prevLoc = 0;
        for (size_t i = 0; i < numElements; ++i)
        {
            AlignmentElementOperation& theOp = theFragment.getAlignmentType(i);
            SequenceLocation& qLoc = theFragment.getQueryLocation(i);
            SequenceLocation& tLoc = theFragment.getTargetLocation(i);
            // Do not need to adjust by theOrigin since subtracting by similar.
            position_type& qLocPos = qLoc.position();
            position_type& tLocPos = tLoc.position();

            // First check for skipped bases on reference.
#ifdef PROG_DEBUG
            assert(tLocPos.lower() >= (prevLoc + 1));
#endif
            if ((i > 0) && (tLocPos.lower() > (prevLoc + 1)))
            {
                long distSkipped = tLocPos.lower() - prevLoc - 1;
                cigarString += convertToString<long>(distSkipped) + "N";
            }
            prevLoc = tLocPos.upper();

            long opSize;
            switch (theOp)
            {
                case MatchOrMismatch :
                case NoChange :
                case SNPs :
                    opSize = tLocPos.upper() - tLocPos.lower() + 1;
                    cigarString += convertToString<long>(opSize);
                    if (useNewCigar)
                    {
                        if (theOp == NoChange)
                        {
                            cigarString += "=";
                        }
                        else if (theOp == SNPs)
                        {
                            cigarString += "X";
                        }
                        else  // MatchOrMismatch
                        {
                            cigarString += "M";
                        }
                    }
                    else
                    {
                        cigarString += "M";    // Use this for all three.
                    }
                    break;
                case Insertion :
                    opSize = qLocPos.upper() - qLocPos.lower() + 1;
                    cigarString += convertToString<long>(opSize) + "I";
                    break;
                case Deletion :
                    opSize = tLocPos.upper() - tLocPos.lower() + 1;
                    cigarString += convertToString<long>(opSize) + "D";
                    break;
                case Padding :
                case SoftClip :
                case HardClip :
                    if (tLoc.info().hasKey("element_silent_length"))
                    {
                        cigarString += tLoc.info().getQualifier("element_silent_length");
                    }
                    else
                    {
                        cigarString += "0";
                    }
                    if (theOp == Padding)
                    {
                        cigarString += "P";
                    }
                    else if (theOp == SoftClip)
                    {
                        cigarString += "S";
                    }
                    else if (theOp == HardClip)
                    {
                        cigarString += "H";
                    }
                    break;
                default :
#ifdef PROG_DEBUG
                    // Should never get here
                    assert(true == false);
#endif
                    break;
            }
        }

        out << cigarString << "\t";

        if (targetInfo.hasKey("MRNM"))
        {
            out << targetInfo.getQualifier("MRNM") << "\t";
        }
        else
        {
            out << "*\t";
        }

        if (targetInfo.hasKey("MPOS"))
        {
            string tempMpos = targetInfo.getQualifier("MPOS");
            long lMpos = convertFromString<long>(tempMpos);
            // Below +1 is because SAM uses 1-based index.
            lMpos += 1 - theOrigin;
            out << lMpos << "\t";
        }
        else
        {
            out << "0\t";
        }

        if (targetInfo.hasKey("ISIZE"))
        {
            out << targetInfo.getQualifier("ISIZE") << "\t";
        }
        else
        {
            out << "0\t";
        }


        size_t theIndex;
        if (printFasta)
        {
            theIndex = ptrReads->getIndex(globalQueryLocation.contig());
            out << ptrReads->getSequence(theIndex) << "\t";
        }
        else
        {
            out << "*\t";
        }
        if (printQuals)
        {
#ifdef PROG_DEBUG
            assert(printFasta);
#endif
            theIndex = ptrReads->getIndex(globalQueryLocation.contig());
            out << ptrQualities->getQual(theIndex);
        }
        else
        {
            out << "*";
        }

        if (targetInfo.hasKey("optional_fields"))
        {
            out << "\t" << targetInfo.getQualifier("optional_fields");
        }

        out << endl;
    }


private:
    ostream& out;
    Fasta *ptrReads;
    Qual *ptrQualities;
    bool printFasta, printQuals, useNewCigar;


};

#endif

