#ifndef MOLBIOLIB_BROADQLTOUTALIGNMENTFRAGMENTWRITER_H
#define MOLBIOLIB_BROADQLTOUTALIGNMENTFRAGMENTWRITER_H 1

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
#include "src/Functions/ReaderWriters/Alignments/BroadQltoutAlignmentFragmentReader.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"




/** \file BroadQltoutAlignmentFragmentWriter.hpp
 * Writes a QUERY line in a presumably open Broad Qltout file given a fragment.
 * The usage is
 * \code
 * SequenceAlignmentFragment theFragment;   // Initialize later...
 * // The fragment positions must be 1-based.
 * // The fragment positions must also be closed-intervals - they will be
 * //    converted to half-open-intervals by the writer.
 * ostream out;   // Initialize later.  Can also be of type Ofstream, etc.
 * BroadQltoutAlignmentFragmentWriter theWriter(out, BroadContigNumber, nullStringQualifiers, BroadContigNumber, nullStringQualifiers);
 * \endcode
 * where the last four parameters a defaults.  One passes the same last four
 * parameters as one would have done for
 * the #BroadQltoutAlignmentFragmentReader.
 *
 * One then can do
 * \code
 * theWriter.write(theFragment);
 * theWriter.writeWithForcedParams(theFragment, false, 0, false, 0);
 * \endcode
 * where the last four parameters are optional.  If specify <code>true</code>,
 * this forces the query or target (respectively) to be written as the number
 * passed.  The default values are shown.
 *
 */
class BroadQltoutAlignmentFragmentWriter : public SequenceAlignmentFragmentWriter
{
public:

    // Can initialize an instance of this class using just the
    // BroadContigTranslation parameters and the ostream.  Then,
    // use the below two to make the reverse qualifiers.  This way,
    // one can make the mapping separately for memory conservation.
    void generateQueryReverseQualifiers(StringQualifiers& theQueryQualifiers,
                                        bool clearOrgQuals = false)
    {
        queryReversedQualifiers.clear();
        theQueryQualifiers.generateReverseQualifiers(queryReversedQualifiers);
        if (clearOrgQuals)
        {
            theQueryQualifiers.clear();
        }
    }
    void generateTargetReverseQualifiers(StringQualifiers& theTargetQualifiers,
                                         bool clearOrgQuals = false)
    {
        targetReversedQualifiers.clear();
        theTargetQualifiers.generateReverseQualifiers(targetReversedQualifiers);
        if (clearOrgQuals)
        {
            theTargetQualifiers.clear();
        }
    }

    // clearOrgQuals parameter included in situations where memory is tight.
    // Will clear the qualifier after generating the reverse.
    void setTrans(BroadContigTranslation inputWhichQueryTrans = BroadContigNumber,
                  StringQualifiers& theQueryQualifiers = nullStringQualifiers,
                  BroadContigTranslation inputWhichTargetTrans = BroadContigNumber,
                  StringQualifiers& theTargetQualifiers = nullStringQualifiers,
                  bool clearOrgQuals = false)
    {
        whichQueryTrans = inputWhichQueryTrans;
        whichTargetTrans = inputWhichTargetTrans;

        if (whichQueryTrans == BroadContigStringQualifiers)
        {
            generateQueryReverseQualifiers(theQueryQualifiers, clearOrgQuals);
        }
        if (whichTargetTrans == BroadContigStringQualifiers)
        {
            generateTargetReverseQualifiers(theTargetQualifiers);
        }
    }


    // Below can be used to construct this class with only the ostream.  Can
    // set up translations later with setTrans(...).
    BroadQltoutAlignmentFragmentWriter(ostream& ofp,
                                       BroadContigTranslation inputWhichQueryTrans = BroadContigNumber,
                                       StringQualifiers& inputTheQueryQualifiers = nullStringQualifiers,
                                       BroadContigTranslation inputWhichTargetTrans = BroadContigNumber,
                                       StringQualifiers& inputTheTargetQualifiers = nullStringQualifiers,
                                       bool clearOrgQuals = false) : out(ofp)
    {
        setTrans(inputWhichQueryTrans, inputTheQueryQualifiers, inputWhichTargetTrans, inputTheTargetQualifiers, clearOrgQuals);
    }


    BroadQltoutAlignmentFragmentWriter(ostream& ofp,
                                       BroadContigTranslation inputWhichQueryTrans = BroadContigNumber,
                                       BroadContigTranslation inputWhichTargetTrans = BroadContigNumber) : out(ofp)
    {
        whichQueryTrans = inputWhichQueryTrans;
        whichTargetTrans = inputWhichTargetTrans;
    }





    void writeWithForcedParams(SequenceAlignmentFragment& theFragment,
                               bool forceQueryNum = false, size_t queryNum = 0,
                               bool forceTargetNum = false, size_t targetNum = 0)
    {
        string defaultChrString = "";
        if ((whichQueryTrans == BroadContigNumberDefaultChr) ||
                (whichTargetTrans == BroadContigNumberDefaultChr))
        {
            defaultChrString = DEFAULT_CONTIG_STRING;
        }
        size_t defaultChrStringLength = defaultChrString.size();
        size_t chrPos;
        size_t broadChrNum;

        out << "QUERY\t";

        SequenceLocation& globalQueryLocation = theFragment.globalQueryLocation;
        SequenceLocation::element_type theOrigin = globalQueryLocation.origin();
        string queryContig = globalQueryLocation.contig();
        string qltoutQuery = "";
        if (!forceQueryNum)
        {
            switch (whichQueryTrans)
            {
                case BroadContigNumber :
                    qltoutQuery = queryContig;
                    break;
                case BroadContigNumberDefaultChr :
                    chrPos = queryContig.find(defaultChrString) + defaultChrStringLength;
                    qltoutQuery = queryContig.substr(chrPos);   // The string after the
                    // default chr string.
                    broadChrNum = convertFromString<size_t>(qltoutQuery);
                    qltoutQuery = convertToString<size_t>(broadChrNum - 1);
                    break;
                case BroadContigStringQualifiers :
#ifdef DEBUG
                    // Must have unique mapping back.
                    assert(queryReversedQualifiers[queryContig].size() == 1);
#endif
                    qltoutQuery = queryReversedQualifiers[queryContig][0];
                    break;
                default :
#ifdef PROG_DEBUG
                    // This should never happen.
                    assert(true == false);
#endif
                    break;
            }
        }
        else
        {
            qltoutQuery = convertToString<size_t>(queryNum);
        }
        out << qltoutQuery << "\t";

        SequenceLocation::position_type queryPos = globalQueryLocation.position();
        StringQualifiers& queryInfo = globalQueryLocation.info();
        // Upper is added by 1 since Broad uses open-ended coordinates.
        out << queryPos.lower() - theOrigin << "\t" << queryPos.upper() - theOrigin + 1 << "\t";
        if (queryInfo.hasKey("query_contig_length"))
        {
            out << queryInfo.getQualifier("query_contig_length") << "\t";    // Last is length.
        }
        else
        {
            out << "0\t";    // No query length found.
        }

        SequenceLocation& globalTargetLocation = theFragment.globalTargetLocation;
        if (globalTargetLocation.strand() == "+")
        {
            out << "0\t";
        }
        else
        {
            out << "1\t";
        }
        string targetContig = globalTargetLocation.contig();
        string qltoutTarget = "";
        if (!forceTargetNum)
        {
            switch (whichTargetTrans)
            {
                case BroadContigNumber :
                    qltoutTarget = targetContig;
                    break;
                case BroadContigNumberDefaultChr :
                    chrPos = targetContig.find(defaultChrString) + defaultChrStringLength;
                    qltoutTarget = targetContig.substr(chrPos);   // The string after the
                    // default chr string.
                    broadChrNum = convertFromString<size_t>(qltoutTarget);
                    qltoutTarget = convertToString<size_t>(broadChrNum - 1);
                    break;
                case BroadContigStringQualifiers :
#ifdef DEBUG
                    // Must have unique mapping back.
                    assert(targetReversedQualifiers[targetContig].size() == 1);
#endif
                    qltoutTarget = targetReversedQualifiers[targetContig][0];
                    break;
                default :
#ifdef PROG_DEBUG
                    // This should never happen.
                    assert(true == false);
#endif
                    break;
            }
        }
        else
        {
            qltoutTarget = convertToString<size_t>(targetNum);
        }
        out << qltoutTarget << "\t";

        SequenceLocation::position_type targetPos = globalTargetLocation.position();
        StringQualifiers& targetInfo = globalTargetLocation.info();
        // Upper is added by 1 since Broad uses open-ended coordinates.
        out << targetPos.lower() - theOrigin << "\t" << targetPos.upper() - theOrigin + 1 << "\t";
        if (targetInfo.hasKey("target_contig_length"))
        {
            out << targetInfo.getQualifier("target_contig_length") << "\t";    // Last is length.
        }
        else
        {
            out << "0\t";    // No target contig length found.
        }


        size_t numFragments = theFragment.size();
        out << numFragments;   // Number of blocks.  No tabs.  Do below.

        SequenceAlignmentFragment::row_type theRow;
        for (size_t i = 0; i < numFragments; ++i)
        {
            AlignmentElementOperation currOp = theFragment.getAlignmentType(i);
            SequenceLocation targetLocation = theFragment.getTargetLocation(i);
            // Do not need to subtract theOrigin since
            // doing subtraction of like below.
            SequenceLocation::position_type targetLocPos = targetLocation.position();
            StringQualifiers& targetInfo = targetLocation.info();
            long absG = 0;
            if (targetInfo.hasKey("element_gapLength"))
            {
                absG = convertFromString<long>(targetInfo.getQualifier("element_gapLength"));
            }
            long e = 0;
            if (targetInfo.hasKey("element_numErrors"))
            {
                e = convertFromString<long>(targetInfo.getQualifier("element_numErrors"));
            }
            long g = absG;
            if (currOp == Insertion)
            {
                g = -absG;
            }
            long b = targetLocPos.upper() - targetLocPos.lower() + 1;
            // Below is Broad-specific.
            // The gap-length does not count towards the block length.
            if (g != 0)
            {
                b = max(static_cast<long>(0), b-absG);

            }
            out << "\t" << g << "\t" << b << "\t" << e;
        }
        out << endl;
    }


    void write(SequenceAlignmentFragment& theFragment)
    {
        writeWithForcedParams(theFragment);
    }


private:
    ostream& out;
    BroadContigTranslation whichQueryTrans, whichTargetTrans;
    StringQualifiers::reverseQualifier_type queryReversedQualifiers, targetReversedQualifiers;

};


#endif
